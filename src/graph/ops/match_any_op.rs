use block_aligner::{cigar::*, scan_block::*, scores::*};


use memchr::memmem;

use thread_local::*;

use std::cell::RefCell;
use std::marker::Send;

use crate::graph::*;
use crate::seed_search::*;
use crate::Patterns;

thread_local! {
    static MATCH_ANY_FLAGS: RefCell<Vec<Option<Option<isize>>>> = RefCell::new(Vec::new());
    static MATCH_ANY_USED: RefCell<Vec<usize>> = RefCell::new(Vec::new());
    static MATCH_ANY_CANDS: RefCell<Vec<(usize, Option<isize>)>> = RefCell::new(Vec::new());
}

pub struct MatchAnyOp {
    required_names: Vec<LabelOrAttr>,
    label: Label,
    new_labels: [Option<Label>; 3],
    patterns: Patterns,
    // Length of literal patterns by index (None for expr-based patterns)
    literal_len: Vec<Option<usize>>,
    max_literal_len: usize,
    all_literals: bool,
    match_type: MatchType,
    aligner: ThreadLocal<Option<RefCell<Box<dyn Aligner + Send>>>>,
    seed_searcher: Option<SeedSearchers>,
}

impl MatchAnyOp {
    const NAME: &'static str = "MatchAnyOp";

    /// Match any one of multiple patterns in an interval.
    ///
    /// Patterns can be arbitrary expressions, so you can use any existing labeled intervals or
    /// attributes as patterns.
    ///
    /// You can also include arbitrary extra attributes for each pattern. The corresponding attributes
    /// for the matched pattern will be stored into the input labeled interval.
    ///
    /// The transform expression must have one input label and the number of output labels is
    /// determined by the [`MatchType`].
    ///
    /// Example `transform_expr` for local-alignment-based pattern matching:
    /// `tr!(seq1.* -> seq1.before, seq1.aligned, seq1.after)`.
    /// The input labeled interval will get a new attribute (`seq1.*.my_patterns`) that is set to the pattern
    /// that is matched. If no pattern matches, then it will be set to false.
    pub fn new(transform_expr: TransformExpr, patterns: Patterns, match_type: MatchType) -> Self {
        let mut new_labels = [None, None, None];

        transform_expr.check_size(1, match_type.num_mappings(), Self::NAME);
        for i in 0..match_type.num_mappings() {
            new_labels[i] = transform_expr.after_label(i, Self::NAME);
        }
        transform_expr.check_same_str_type(Self::NAME);

        let seed_searcher = Self::get_searcher(&patterns, &match_type);
        let max_literal_len = patterns
            .iter_literals()
            .map(|(_, p)| p.len())
            .max()
            .unwrap_or(0);
        let all_literals = patterns.iter_exprs().count() == 0;
        // Build a map of literal lengths by pattern index to allow quick sorting of candidates
        let mut literal_len = vec![None; patterns.patterns().len()];
        for (i, p) in patterns.iter_literals() {
            if i < literal_len.len() {
                literal_len[i] = Some(p.len());
            }
        }
        let mut required_names = vec![transform_expr.before(0).into()];
        required_names.extend(
            patterns
                .iter_exprs()
                .flat_map(|(_, e)| e.required_names().into_iter()),
        );

        Self {
            required_names,
            label: transform_expr.before(0),
            new_labels,
            patterns,
            literal_len,
            max_literal_len,
            all_literals,
            match_type,
            aligner: ThreadLocal::new(),
            seed_searcher,
        }
    }

    fn get_searcher(patterns: &Patterns, match_type: &MatchType) -> Option<SeedSearchers> {
        let min_len = patterns
            .iter_literals()
            .map(|(_, p)| p.len())
            .min()
            .unwrap_or(0);
        let k = match_type.k(min_len);

        use SeedSearchers::*;
        let res = match k {
            0..=1 => return None,
            2 => SmallSearcher::<2>::new(patterns.iter_literals()).map(Small2),
            3 => SmallSearcher::<3>::new(patterns.iter_literals()).map(Small3),
            4 => SmallSearcher::<4>::new(patterns.iter_literals()).map(Small4),
            5 => SmallSearcher::<5>::new(patterns.iter_literals()).map(Small5),
            6 => SmallSearcher::<6>::new(patterns.iter_literals()).map(Small6),
            _ => Err(()),
        };

        if let Ok(s) = res {
            Some(s)
        } else {
            Some(General(GeneralSearcher::new(patterns.iter_literals(), k)))
        }
    }
}

impl<T: crate::trace::Trace> GraphNode<T> for MatchAnyOp {
    fn run_inner(&self, mut read: Read) -> Result<(Option<Read>, bool)> {
        let text = read
            .substring(self.label.str_type, self.label.label)
            .map_err(|e| Error::NameError {
                source: e,
                read: read.clone(),
                context: Self::NAME,
            })?;

        use MatchType::*;
        let aligner_cell = self.aligner.get_or(|| {
            let init_len = if self.max_literal_len > 0 {
                self.max_literal_len * 2
            } else {
                text.len() * 2
            };

            match self.match_type {
                GlobalAln(_) => Some(RefCell::new(Box::new(GlobalLocalAligner::<false>::new(
                    init_len,
                )))),
                LocalAln { .. } => Some(RefCell::new(Box::new(GlobalLocalAligner::<true>::new(
                    init_len,
                )))),
                PrefixAln { .. } => Some(RefCell::new(Box::new(PrefixSuffixAligner::<true>::new(
                    init_len,
                )))),
                SuffixAln { .. } => Some(RefCell::new(Box::new(
                    PrefixSuffixAligner::<false>::new(init_len),
                ))),
                _ => None,
            }
        });

        let additional = |identity: f64, pattern_len: usize| {
            ((1.0 - identity).max(0.0) * (pattern_len as f64)).ceil() as usize
        };

        // Reuse thread-local buffers to record candidate pattern indices and their optional offsets.
        // flags[i] = Some(Some(text_i)) means candidate with offset; Some(None) means candidate with no offset; None means absent.
        let patterns_len = self.patterns.patterns().len();

        if let Some(seed_searcher) = &self.seed_searcher {
            let (text_slice, text_offset, use_i) = match self.match_type {
                Exact => (text, 0, false),
                ExactPrefix => (&text[..text.len().min(self.max_literal_len)], 0, false),
                ExactSuffix => {
                    let offset = text.len().saturating_sub(self.max_literal_len);
                    (&text[offset..], offset, false)
                }
                ExactSearch => (text, 0, true),
                ExactBoundedMatch { from, to } => {
                    let to = text.len().min(to);
                    (&text[from..to], 0, false)
                },
                Hamming(_) => (text, 0, false),
                HammingPrefix(_) => (&text[..text.len().min(self.max_literal_len)], 0, false),
                HammingSuffix(_) => {
                    let offset = text.len().saturating_sub(self.max_literal_len);
                    (&text[offset..], offset, false)
                }
                HammingSearch(_) => (text, 0, true),
                HammingBoundedMatch {
                    threshold: _,
                    from,
                    to,
                } => {
                    let to = text.len().min(to);
                    (&text[from..to], 0, false)
                },
                GlobalAln(_) => (text, 0, false),
                LocalAln { .. } => (text, 0, true),
                PrefixAln { identity, .. } => (
                    &text[..text
                        .len()
                        .min(self.max_literal_len + additional(identity, self.max_literal_len))],
                    0,
                    false,
                ),
                SuffixAln { identity, .. } => {
                    let offset = text.len().saturating_sub(
                        self.max_literal_len + additional(identity, self.max_literal_len),
                    );
                    (&text[offset..], offset, false)
                }
            };

            MATCH_ANY_FLAGS.with(|flags_cell| {
                MATCH_ANY_USED.with(|used_cell| {
                    let mut flags_ref = flags_cell.borrow_mut();
                    let mut used_ref = used_cell.borrow_mut();
                    if flags_ref.len() < patterns_len { flags_ref.resize(patterns_len, None); }
                    seed_searcher.search(text_slice, |SeedMatch { pattern_idx, pattern_i, text_i }| {
                        let ti = if use_i {
                            Some(((text_offset + text_i) as isize) - (pattern_i as isize))
                        } else {
                            None
                        };
                        if flags_ref[pattern_idx].is_none() { used_ref.push(pattern_idx); }
                        flags_ref[pattern_idx] = Some(ti);
                    });
                });
            });
        } else {
            MATCH_ANY_FLAGS.with(|flags_cell| {
                MATCH_ANY_USED.with(|used_cell| {
                    let mut flags_ref = flags_cell.borrow_mut();
                    let mut used_ref = used_cell.borrow_mut();
                    if flags_ref.len() < patterns_len { flags_ref.resize(patterns_len, None); }
                    used_ref.clear();
                    for (i, _) in self.patterns.iter_literals() {
                        if flags_ref[i].is_none() { used_ref.push(i); }
                        flags_ref[i] = Some(None);
                    }
                });
            });
        }

        if !self.all_literals {
            MATCH_ANY_FLAGS.with(|flags_cell| {
                MATCH_ANY_USED.with(|used_cell| {
                    let mut flags_ref = flags_cell.borrow_mut();
                    let mut used_ref = used_cell.borrow_mut();
                    if flags_ref.len() < patterns_len { flags_ref.resize(patterns_len, None); }
                    for (i, _) in self.patterns.iter_exprs() {
                        if flags_ref[i].is_none() { used_ref.push(i); }
                        flags_ref[i] = Some(None);
                    }
                });
            });
        }

        // Collect to a TLS vector of (idx, text_i), then move it out to avoid per-read allocations.
        let mut seed_hits_vec: Vec<(usize, Option<isize>)> = MATCH_ANY_CANDS.with(|cands_cell| {
            MATCH_ANY_FLAGS.with(|flags_cell| {
                MATCH_ANY_USED.with(|used_cell| {
                    let mut cands = cands_cell.borrow_mut();
                    cands.clear();
                    let mut flags_ref = flags_cell.borrow_mut();
                    let mut used_ref = used_cell.borrow_mut();
                    cands.reserve(used_ref.len());
                    for &idx in used_ref.iter() {
                        let ti = flags_ref[idx].take().unwrap();
                        cands.push((idx, ti));
                    }
                    used_ref.clear();
                    std::mem::take(&mut *cands)
                })
            })
        });
        // Heuristic: prioritize longer literal patterns to increase early pruning by max_matches.
        if seed_hits_vec.len() > 16 {
            seed_hits_vec.sort_by_key(|(idx, _)| std::cmp::Reverse(self.literal_len.get(*idx).and_then(|x| *x).unwrap_or(0)));
        }

        let mut max_matches = 0;
        let mut max_pattern = None;
        let mut max_pattern_idx = std::usize::MAX;
        let mut max_cut_pos1 = 0;
        let mut max_cut_pos2 = 0;
        let mut multimatches = false;

        'candidates: for (pattern_idx, text_i) in seed_hits_vec.iter().copied() {
            let pattern = &self.patterns.patterns()[pattern_idx];
            let pattern_str_cow = pattern.get(&read).map_err(|e| Error::NameError {
                source: e,
                read: read.clone(),
                context: Self::NAME,
            })?;
            let pattern_str: &[u8] = &pattern_str_cow;
            let pattern_len = pattern_str.len();

            if max_matches > pattern_len {
                continue;
            }

            let matches = match self.match_type {
                Exact => {
                    if text == pattern_str {
                        Some((pattern_len, pattern_len, 0))
                    } else {
                        None
                    }
                }
                ExactPrefix => {
                    if pattern_len <= text.len() && &text[..pattern_len] == pattern_str {
                        Some((pattern_len, pattern_len, 0))
                    } else {
                        None
                    }
                }
                ExactSuffix => {
                    if pattern_len <= text.len() && &text[text.len() - pattern_len..] == pattern_str
                    {
                        Some((pattern_len, text.len() - pattern_len, 0))
                    } else {
                        None
                    }
                }
                ExactSearch => {
                    let (text_start, text_end) = if let Some(text_i) = text_i {
                        (
                            text_i.max(0) as usize,
                            text.len().min((text_i + (pattern_len as isize)) as usize),
                        )
                    } else {
                        (0, text.len())
                    };
                    let text_around = &text[text_start..text_end];
                    memmem::find(text_around, pattern_str)
                        .map(|i| (pattern_len, text_start + i, text_start + i + pattern_len))
                }
                ExactBoundedMatch { from, to } => {
                    let to = text.len().min(to);
                    let text_around = &text[from..to];
                    memmem::find(text_around, pattern_str)
                        .map(|i| (pattern_len, from + i, from + i + pattern_len))
                }
                Hamming(t) => {
                    let t = t.get(pattern_len);
                    hamming(text, pattern_str, t).map(|m| (m, pattern_len, 0))
                }
                HammingPrefix(t) => {
                    if pattern_len <= text.len() {
                        let t = t.get(pattern_len);
                        hamming(&text[..pattern_len], pattern_str, t).map(|m| (m, pattern_len, 0))
                    } else {
                        None
                    }
                }
                HammingSuffix(t) => {
                    if pattern_len <= text.len() {
                        let t = t.get(pattern_len);
                        hamming(&text[text.len() - pattern_len..], pattern_str, t)
                            .map(|m| (m, text.len() - pattern_len, 0))
                    } else {
                        None
                    }
                }
                HammingSearch(t) => {
                    let (text_start, text_end) = if let Some(text_i) = text_i {
                        (
                            text_i.max(0) as usize,
                            text.len().min((text_i + (pattern_len as isize)) as usize),
                        )
                    } else {
                        (0, text.len())
                    };
                    let text_around = &text[text_start..text_end];
                    let t = t.get(pattern_len);
                    hamming_search(text_around, pattern_str, t).map(|(m, start_idx, end_idx)| {
                        (m, text_start + start_idx, text_start + end_idx)
                    })
                }
                HammingBoundedMatch {
                    threshold: t,
                    from,
                    to,
                } => {
                    let t = t.get(pattern_len);
                    let to = text.len().min(to);
                    let text_around = &text[from..to];
                    hamming_search(text_around, pattern_str, t)
                        .map(|(m, start_idx, end_idx)| (m, from + start_idx, from + end_idx))
                }
                GlobalAln(identity) => {
                    if (identity - 1.0).abs() < f64::EPSILON {
                        if text == pattern_str {
                            Some((pattern_len, pattern_len, 0))
                        } else {
                            None
                        }
                    } else if text.len() == pattern_len {
                        // Fast path: equal lengths. Use Hamming to accept early if identity threshold is met
                        // without invoking the aligner (gapless alignment is valid for global).
                        let thr = (identity * (pattern_len as f64)).ceil() as usize;
                        hamming(text, pattern_str, thr).map(|m| (m, pattern_len, 0))
                    } else {
                        aligner_cell
                            .as_ref()
                            .unwrap()
                            .borrow_mut()
                            .align(text, pattern_str, identity, identity)
                            .map(|(m, _, end_idx)| (m, end_idx, 0))
                    }
                }
                LocalAln { identity, overlap } => {
                    if (identity - 1.0).abs() < f64::EPSILON && (overlap - 1.0).abs() < f64::EPSILON {
                        if let Some(start) = memmem::find(text, pattern_str) {
                            Some((pattern_len, start, start + pattern_len))
                        } else {
                            None
                        }
                    } else {
                        let a = additional(identity, pattern_len) as isize;
                        let (text_start, text_end) = if let Some(text_i) = text_i {
                            (
                                (text_i - a).max(0) as usize,
                                text.len().min((text_i + (pattern_len as isize) + a) as usize),
                            )
                        } else {
                            (0, text.len())
                        };
                        let text_around = &text[text_start..text_end];
                        aligner_cell
                            .as_ref()
                            .unwrap()
                            .borrow_mut()
                            .align(text_around, pattern_str, identity, overlap)
                            .map(|(m, start_idx, end_idx)| (m, text_start + start_idx, text_start + end_idx))
                    }
                }
                PrefixAln { identity, overlap } => {
                    if (identity - 1.0).abs() < f64::EPSILON && (overlap - 1.0).abs() < f64::EPSILON {
                        if pattern_len <= text.len() && &text[..pattern_len] == pattern_str {
                            Some((pattern_len, pattern_len, 0))
                        } else {
                            None
                        }
                    } else {
                        let a = additional(identity, pattern_len);
                        aligner_cell
                            .as_ref()
                            .unwrap()
                            .borrow_mut()
                            .align(&text[..text.len().min(pattern_len + a)], pattern_str, identity, overlap)
                            .map(|(m, _, end_idx)| (m, end_idx, 0))
                    }
                }
                SuffixAln { identity, overlap } => {
                    if (identity - 1.0).abs() < f64::EPSILON && (overlap - 1.0).abs() < f64::EPSILON {
                        if pattern_len <= text.len() && &text[text.len() - pattern_len..] == pattern_str {
                            Some((pattern_len, text.len() - pattern_len, 0))
                        } else {
                            None
                        }
                    } else {
                        let a = additional(identity, pattern_len);
                        let text_start = text.len().saturating_sub(pattern_len + a);
                        aligner_cell
                            .as_ref()
                            .unwrap()
                            .borrow_mut()
                            .align(&text[text_start..], pattern_str, identity, overlap)
                            .map(|(m, start_idx, _)| (m, text_start + start_idx, 0))
                    }
                }
            };

            if let Some((matches, cut_pos1, cut_pos2)) = matches {
                if matches > max_matches {
                    max_matches = matches;
                    max_pattern = Some((pattern_str_cow, pattern.attrs()));
                    max_pattern_idx = pattern_idx;
                    max_cut_pos1 = cut_pos1;
                    max_cut_pos2 = cut_pos2;
                    multimatches = false;
                    // Early-exit: if all patterns are literals and we achieved a perfect match
                    // of the maximum literal length, no other candidate can surpass this.
                    if self.all_literals && matches == self.max_literal_len {
                        break 'candidates;
                    }
                } else if matches == max_matches && pattern_idx != max_pattern_idx {
                    multimatches = true;
                }
            }
        }

        // Return the candidates Vec to TLS pool for reuse
        MATCH_ANY_CANDS.with(|cell| {
            *cell.borrow_mut() = seed_hits_vec;
        });

        if let Some((pattern_str_cow, pattern_attrs)) = max_pattern {
            // Convert only if we actually need to store the pattern bytes; ensure no immutable borrow of `read` remains
            // before we take a mutable borrow for mapping updates.
            let maybe_owned = if self.patterns.pattern_name().is_some() {
                Some(pattern_str_cow.into_owned())
            } else {
                // Explicitly drop to end borrow.
                drop(pattern_str_cow);
                None
            };
            let mapping = read
                .mapping_mut(self.label.str_type, self.label.label)
                .unwrap();

            if let Some(pattern_name) = self.patterns.pattern_name() {
                // safe to unwrap: we created owned bytes when pattern_name was Some
                *mapping.data_mut(pattern_name) = Data::Bytes(maybe_owned.unwrap());
            }

            if let Some(multimatch_name) = self.patterns.multimatch_name() {
                *mapping.data_mut(multimatch_name) = Data::Bool(multimatches);
            }

            for (&attr, data) in self.patterns.attr_names().iter().zip(pattern_attrs) {
                *mapping.data_mut(attr) = data.clone();
            }

            match self.match_type.num_mappings() {
                1 => {
                    let start = mapping.start;
                    let str_mappings = read.str_mappings_mut(self.label.str_type).unwrap();
                    str_mappings.add_mapping(
                        self.new_labels[0].as_ref().map(|l| l.label),
                        start,
                        max_cut_pos1,
                    );
                }
                2 => {
                    read.cut(
                        self.label.str_type,
                        self.label.label,
                        self.new_labels[0].as_ref().map(|l| l.label),
                        self.new_labels[1].as_ref().map(|l| l.label),
                        max_cut_pos1 as isize,
                    )
                    .unwrap_or_else(|e| panic!("Error in {}: {e}", Self::NAME));
                }
                3 => {
                    let offset = mapping.start;
                    let mapping_len = mapping.len;

                    let str_mappings = read.str_mappings_mut(self.label.str_type).unwrap();
                    str_mappings.add_mapping(
                        self.new_labels[0].as_ref().map(|l| l.label),
                        offset,
                        max_cut_pos1,
                    );
                    str_mappings.add_mapping(
                        self.new_labels[1].as_ref().map(|l| l.label),
                        offset + max_cut_pos1,
                        max_cut_pos2 - max_cut_pos1,
                    );
                    str_mappings.add_mapping(
                        self.new_labels[2].as_ref().map(|l| l.label),
                        offset + max_cut_pos2,
                        mapping_len - max_cut_pos2,
                    );
                }
                _ => unreachable!(),
            }
        } else {
            if let Some(pattern_name) = self.patterns.pattern_name() {
                *read
                    .mapping_mut(self.label.str_type, self.label.label)
                    .unwrap()
                    .data_mut(pattern_name) = Data::Bool(false);
            }
        }

        Ok((Some(read), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}

fn hamming(a: &[u8], b: &[u8], threshold: usize) -> Option<usize> {
    if a.len() != b.len() {
        return None;
    }

    let a_ptr = a.as_ptr();
    let b_ptr = b.as_ptr();
    let n = a.len();
    let mut res = 0;
    let mut i = 0;

    unsafe {
        while i < (n / 8) * 8 {
            let a_word = std::ptr::read_unaligned(a_ptr.add(i) as *const u64);
            let b_word = std::ptr::read_unaligned(b_ptr.add(i) as *const u64);

            let xor = a_word ^ b_word;
            let or1 = xor | (xor >> 1);
            let or2 = or1 | (or1 >> 2);
            let or3 = or2 | (or2 >> 4);
            let mask = or3 & 0x0101010101010101u64;
            res += mask.count_ones() as usize;

            i += 8;
        }

        if i < n {
            let a_word = read_rest_u64(a_ptr.add(i), n - i);
            let b_word = read_rest_u64(b_ptr.add(i), n - i);

            let xor = a_word ^ b_word;
            let or1 = xor | (xor >> 1);
            let or2 = or1 | (or1 >> 2);
            let or3 = or2 | (or2 >> 4);
            let mask = or3 & 0x0101010101010101u64;
            res += mask.count_ones() as usize;
        }
    }

    let matches = n - res;

    if matches >= threshold {
        Some(matches)
    } else {
        None
    }
}

unsafe fn read_rest_u64(ptr: *const u8, len: usize) -> u64 {
    let addr = ptr as usize;
    let start_page = addr >> 12;
    let end_page = (addr + 7) >> 12;

    if start_page == end_page {
        std::ptr::read_unaligned(ptr as *const u64) & ((1u64 << (len * 8)) - 1)
    } else {
        let mut res = 0u64;
        let mut i = 0;

        while i < len {
            res |= (*ptr.add(i) as u64) << (i * 8);
            i += 1;
        }

        res
    }
}

fn hamming_search(a: &[u8], b: &[u8], threshold: usize) -> Option<(usize, usize, usize)> {
    let mut best_match = None;

    for (i, w) in a.windows(b.len()).enumerate() {
        if let Some(matches) = hamming(w, b, threshold) {
            if let Some((best_matches, _, _)) = best_match {
                if matches <= best_matches {
                    continue;
                }
            }

            best_match = Some((matches, i, i + b.len()));
        }
    }

    best_match
}

trait Aligner {
    fn align(
        &mut self,
        read: &[u8],
        pattern: &[u8],
        identity_threshold: f64,
        overlap_threshold: f64,
    ) -> Option<(usize, usize, usize)>;
}

struct GlobalLocalAligner<const LOCAL: bool> {
    read_padded: PaddedBytes,
    pattern_padded: PaddedBytes,
    matrix: NucMatrix,
    // always store trace
    block: Block<true, LOCAL, LOCAL, false>,
    cigar: Cigar,
    len: usize,
}

impl<const LOCAL: bool> GlobalLocalAligner<LOCAL> {
    const MAX_SIZE: usize = 512;
    const GAPS: Gaps = Gaps {
        open: -2,
        extend: -1,
    };

    pub fn new(len: usize) -> Self {
        let read_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
        let pattern_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
        let matrix = NucMatrix::new_simple(1, -1);

        let block = Block::<true, LOCAL, LOCAL, false>::new(len, len, Self::MAX_SIZE);
        let cigar = Cigar::new(len, len);

        Self {
            read_padded,
            pattern_padded,
            matrix,
            block,
            cigar,
            len,
        }
    }

    fn resize_if_needed(&mut self, len: usize) {
        if len > self.len {
            self.read_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
            self.pattern_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
            self.block = Block::<true, LOCAL, LOCAL, false>::new(len, len, Self::MAX_SIZE);
            self.cigar = Cigar::new(len, len);
            self.len = len;
        }
    }
}

unsafe impl<const LOCAL: bool> Send for GlobalLocalAligner<LOCAL> {}

impl<const LOCAL: bool> Aligner for GlobalLocalAligner<LOCAL> {
    fn align(
        &mut self,
        read: &[u8],
        pattern: &[u8],
        identity_threshold: f64,
        overlap_threshold: f64,
    ) -> Option<(usize, usize, usize)> {
        self.resize_if_needed(pattern.len().max(read.len()));

        let max_size = pattern
            .len()
            .min(read.len())
            .min(Self::MAX_SIZE);

        self.read_padded.set_bytes::<NucMatrix>(read, max_size);
        self.pattern_padded
            .set_bytes::<NucMatrix>(pattern, max_size);

        let min_size = max_size;

        self.block.align(
            &self.pattern_padded,
            &self.read_padded,
            &self.matrix,
            Self::GAPS,
            min_size..=max_size,
            pattern.len() as i32,
        );

        let res = self.block.res();
        self.block.trace().cigar_eq(
            &self.pattern_padded,
            &self.read_padded,
            res.query_idx,
            res.reference_idx,
            &mut self.cigar,
        );

        let mut matches = 0;
        let mut total = 0;

        // Compute start index only when LOCAL == true; for global alignment (LOCAL == false)
        // the caller ignores the start index, so avoid reverse + index math.
        let mut read_start_idx = if LOCAL { res.reference_idx } else { 0 };
        if LOCAL {
            self.cigar.reverse();
            for i in 0..self.cigar.len() {
                let OpLen { op, len } = self.cigar.get(i);
                match op {
                    Operation::Eq => {
                        read_start_idx -= len;
                        matches += len;
                    }
                    Operation::X | Operation::D => {
                        read_start_idx -= len;
                    }
                    _ => (),
                }
                total += len;
            }
        } else {
            // LOCAL == false: we still need matches/total for scoring, but no need to reverse
            // or compute the start index.
            for i in 0..self.cigar.len() {
                let OpLen { op, len } = self.cigar.get(i);
                if op == Operation::Eq { matches += len; }
                if let Operation::Eq | Operation::X | Operation::D = op { total += len; }
            }
        }

        let identity = (matches as f64) / (total as f64);
        let overlap = (matches as f64) / (pattern.len() as f64);

        if identity >= identity_threshold && overlap >= overlap_threshold {
            Some((matches, read_start_idx, res.reference_idx))
        } else {
            None
        }
    }
}

struct PrefixSuffixAligner<const PREFIX: bool> {
    read_padded: PaddedBytes,
    pattern_padded: PaddedBytes,
    matrix: NucMatrix,
    // always store trace
    block1: Block<true, true, false, true>,  // X-drop
    block2: Block<true, false, false, true>, // no X-drop
    cigar: Cigar,
    len: usize,
}

impl<const PREFIX: bool> PrefixSuffixAligner<PREFIX> {
    const MAX_SIZE: usize = 512;
    const GAPS: Gaps = Gaps {
        open: -2,
        extend: -1,
    };

    pub fn new(len: usize) -> Self {
        let read_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
        let pattern_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
        let matrix = NucMatrix::new_simple(1, -1);

        let block1 = Block::<true, true, false, true>::new(len, len, Self::MAX_SIZE);
        let block2 = Block::<true, false, false, true>::new(len, len, Self::MAX_SIZE);
        let cigar = Cigar::new(len, len);

        Self {
            read_padded,
            pattern_padded,
            matrix,
            block1,
            block2,
            cigar,
            len,
        }
    }

    fn resize_if_needed(&mut self, len: usize) {
        if len > self.len {
            self.read_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
            self.pattern_padded = PaddedBytes::new::<NucMatrix>(len, Self::MAX_SIZE);
            self.block1 = Block::<true, true, false, true>::new(len, len, Self::MAX_SIZE);
            self.block2 = Block::<true, false, false, true>::new(len, len, Self::MAX_SIZE);
            self.cigar = Cigar::new(len, len);
            self.len = len;
        }
    }
}

unsafe impl<const PREFIX: bool> Send for PrefixSuffixAligner<PREFIX> {}

impl<const PREFIX: bool> Aligner for PrefixSuffixAligner<PREFIX> {
    fn align(
        &mut self,
        read: &[u8],
        pattern: &[u8],
        identity_threshold: f64,
        overlap_threshold: f64,
    ) -> Option<(usize, usize, usize)> {
        self.resize_if_needed(pattern.len().max(read.len()));

        let max_size = pattern
            .len()
            .min(read.len())
            .min(Self::MAX_SIZE);

        if PREFIX {
            // reverse sequences to convert to aligning suffix
            self.read_padded.set_bytes_rev::<NucMatrix>(read, max_size);
            self.pattern_padded
                .set_bytes_rev::<NucMatrix>(pattern, max_size);
        } else {
            self.read_padded.set_bytes::<NucMatrix>(read, max_size);
            self.pattern_padded
                .set_bytes::<NucMatrix>(pattern, max_size);
        }

        // first align to get where the pattern starts in the read
        // note that the start gaps in the pattern are free and the alignment
        // can end whenever due to X-drop
        self.block1.align(
            &self.pattern_padded,
            &self.read_padded,
            &self.matrix,
            Self::GAPS,
            max_size..=max_size,
            pattern.len() as i32,
        );

        let res = self.block1.res();
        self.block1.trace().cigar_eq(
            &self.pattern_padded,
            &self.read_padded,
            res.query_idx,
            res.reference_idx,
            &mut self.cigar,
        );

        // use traceback to compute where the alignment started
        let mut read_start_idx = res.reference_idx;
        for i in 0..self.cigar.len() {
            let OpLen { op, len } = self.cigar.get(i);
            match op {
                Operation::Eq | Operation::X | Operation::D => read_start_idx -= len,
                _ => (),
            }
        }

        // skip second alignment if first alignment reaches the end of the read
        if res.reference_idx < read.len() {
            // get the overlapping prefix/suffix region
            if PREFIX {
                self.read_padded
                    .set_bytes::<NucMatrix>(&read[..read.len() - read_start_idx], max_size);
                self.pattern_padded
                    .set_bytes::<NucMatrix>(pattern, max_size);
            } else {
                self.read_padded
                    .set_bytes_rev::<NucMatrix>(&read[read_start_idx..], max_size);
                self.pattern_padded
                    .set_bytes_rev::<NucMatrix>(pattern, max_size);
            }

            // align again with read and pattern switched and reversed so that end gaps in the read
            // are free and the alignment ends at read_start_idx and spans the entire pattern
            self.block2.align(
                &self.read_padded,
                &self.pattern_padded,
                &self.matrix,
                Self::GAPS,
                max_size..=max_size,
                pattern.len() as i32,
            );

            let res = self.block2.res();
            self.block2.trace().cigar_eq(
                &self.read_padded,
                &self.pattern_padded,
                res.query_idx,
                res.reference_idx,
                &mut self.cigar,
            );
        }

        // count matches and total columns for calculating identity and overlap
        // We must include all alignment columns (Eq/X/ins/del) in total regardless of whether the
        // second alignment spanned the entire pattern.
        let mut matches = 0;
        let mut total = 0;
        for i in 0..self.cigar.len() {
            let OpLen { op, len } = self.cigar.get(i);
            if op == Operation::Eq {
                matches += len;
            }
            total += len;
        }

        let identity = (matches as f64) / (total as f64);
        let overlap = (matches as f64) / (pattern.len() as f64);

        if identity >= identity_threshold && overlap >= overlap_threshold {
            let start_idx = if PREFIX { 0 } else { read_start_idx };
            let end_idx = if PREFIX {
                read.len() - read_start_idx
            } else {
                read.len()
            };

            Some((matches, start_idx, end_idx))
        } else {
            None
        }
    }
}
