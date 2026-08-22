#[cfg(not(feature = "seed-baseline"))]
use hashbrown::HashTable;
#[cfg(feature = "seed-baseline")]
use rustc_hash::FxHashMap;

use cfg_if::cfg_if;

pub trait SeedSearcher {
    fn search(&self, text: &[u8], candidate_fn: impl FnMut(SeedMatch));
}

pub enum SeedSearchers {
    Small2(SmallSearcher<2>),
    Small3(SmallSearcher<3>),
    Small4(SmallSearcher<4>),
    Small5(SmallSearcher<5>),
    Small6(SmallSearcher<6>),
    General(GeneralSearcher),
}

impl SeedSearchers {
    pub fn search(&self, text: &[u8], candidate_fn: impl FnMut(SeedMatch)) {
        use SeedSearchers::*;
        match self {
            Small2(s) => s.search(text, candidate_fn),
            Small3(s) => s.search(text, candidate_fn),
            Small4(s) => s.search(text, candidate_fn),
            Small5(s) => s.search(text, candidate_fn),
            Small6(s) => s.search(text, candidate_fn),
            General(s) => s.search(text, candidate_fn),
        }
    }
}

pub struct SmallSearcher<const K: usize> {
    #[allow(unused)]
    pattern_luts: [Aligned<64>; K],
    #[allow(unused)]
    pattern_idxs: Vec<(u32, u32)>,
}

impl<const K: usize> SmallSearcher<K> {
    pub fn new<'a>(
        #[allow(unused)] patterns: impl Iterator<Item = (usize, &'a [u8])>,
    ) -> Result<Self, ()> {
        cfg_if! {
            if #[cfg(any(target_feature = "avx2", all(target_arch = "x86_64", feature = "release-simd")))] {
                const REPEAT: Aligned<64> = Aligned::<64>([0u8; 64]);
                let mut pattern_luts = [REPEAT; K];
                let mut pattern_idxs = Vec::new();
                let mut idx = 0;

                for (pattern_idx, pattern) in patterns {
                    for (pattern_i, kmer) in pattern.windows(K).enumerate() {
                        if idx >= 8 {
                            return Err(());
                        }

                        for (i, &c) in kmer.iter().enumerate() {
                            let c = if i % 2 == 0 { (c as usize) & 0b0000_1111 } else { ((c as usize) >> 1) & 0b0000_1111 };

                            for j in 0..4 {
                                pattern_luts[i].0[c + j * 16] |= (1 << idx) as u8;
                            }
                        }

                        pattern_idxs.push((pattern_idx as u32, pattern_i as u32));
                        idx += 1;
                    }
                }

                Ok(Self { pattern_luts, pattern_idxs })
            } else {
                Err(())
            }
        }
    }
}

impl<const K: usize> SeedSearcher for SmallSearcher<{ K }> {
    fn search(
        &self,
        #[allow(unused)] text: &[u8],
        #[allow(unused)] mut candidate_fn: impl FnMut(SeedMatch),
    ) {
        cfg_if! {
            if #[cfg(any(target_feature = "avx2", all(target_arch = "x86_64", feature = "release-simd")))] {
                const L: usize = 32;

                #[inline(always)]
                unsafe fn inner<const K: usize, const REST: bool>(pattern_luts: &[Aligned<64>; K], pattern_idxs: &[(u32, u32)], text: &[u8], candidate_fn: &mut impl FnMut(SeedMatch), i: usize, len: usize) {
                    let lo_4_mask = _mm256_set1_epi8(0b0000_1111i8);
                    let mut set = _mm256_set1_epi8(-1i8);

                    let load_mask = if REST {
                        load_mask_avx2(len)
                    } else {
                        _mm256_set1_epi8(-1i8)
                    };

                    for (j, pattern_lut) in pattern_luts.iter().enumerate() {
                        let mut chars = if REST {
                            read_rest_avx2(text.as_ptr().add(i + j), len)
                        } else {
                            _mm256_loadu_si256(text.as_ptr().add(i + j) as _)
                        };

                        if j % 2 != 0 {
                            chars = _mm256_srli_epi16(chars, 1);
                        }

                        let pattern_lut = _mm256_load_si256(pattern_lut.0.as_ptr() as _);
                        let curr_set = _mm256_shuffle_epi8(pattern_lut, _mm256_and_si256(chars, lo_4_mask));
                        set = _mm256_and_si256(set, curr_set);
                    }

                    set = _mm256_and_si256(set, load_mask);
                    // Each byte is an unsigned bit set of matching pattern
                    // seeds. A signed greater-than comparison loses values
                    // with bit 7 set, so derive `!= 0` from equality instead.
                    let zero = _mm256_cmpeq_epi8(set, _mm256_setzero_si256());
                    let mut nonzero_mask = !(_mm256_movemask_epi8(zero) as u32);

                    if nonzero_mask > 0 {
                        let mut a = Aligned::<{ L }>([0u8; L]);
                        _mm256_store_si256(a.0.as_mut_ptr() as _, set);

                        while nonzero_mask > 0 {
                            let idx = nonzero_mask.trailing_zeros() as usize;
                            let mut pattern_mask = *a.0.as_ptr().add(idx);
                            while pattern_mask != 0 {
                                let hash_idx = pattern_mask.trailing_zeros() as usize;
                                let (pattern_idx, pattern_i) = *pattern_idxs.as_ptr().add(hash_idx);
                                candidate_fn(SeedMatch { pattern_idx: pattern_idx as usize, pattern_i: pattern_i as usize, text_i: i + idx });
                                pattern_mask &= pattern_mask - 1;
                            }

                            nonzero_mask &= nonzero_mask - 1;
                        }
                    }
                }

                let mut i = 0;

                while i + (K - 1) + L <= text.len() {
                    unsafe {
                        inner::<{ K }, false>(&self.pattern_luts, &self.pattern_idxs, text, &mut candidate_fn, i, L);
                    }
                    i += L;
                }

                if i + (K - 1) < text.len() {
                    unsafe {
                        inner::<{ K }, true>(&self.pattern_luts, &self.pattern_idxs, text, &mut candidate_fn, i, text.len() - i - (K - 1));
                    }
                }
            }
        }
    }
}

pub struct GeneralSearcher {
    k: usize,
    table: HashToPatternIdx,
}

impl GeneralSearcher {
    const B: usize = 8;

    pub fn new<'a>(patterns: impl Iterator<Item = (usize, &'a [u8])>, k: usize) -> Self {
        let mut hashes = Vec::new();

        for (pattern_idx, p) in patterns {
            unsafe {
                Self::get_hashes(p, k, |curr_hashes, len, pattern_i| {
                    hashes.extend(
                        curr_hashes[..len]
                            .iter()
                            .enumerate()
                            .map(|(i, &h)| (h, pattern_idx, pattern_i + i)),
                    );
                });
            }
        }

        Self {
            k,
            table: HashToPatternIdx::new(hashes),
        }
    }

    #[inline(always)]
    unsafe fn get_hashes(s: &[u8], k: usize, mut f: impl FnMut([u64; Self::B], usize, usize)) {
        if s.len() < k {
            return;
        }

        let mut hash = 0u64;
        let mut i = 0;
        let ptr = s.as_ptr();

        while i < k {
            hash = hash.rotate_left(1) ^ wyhash_byte(*ptr.add(i));
            i += 1;
        }

        while i + Self::B <= s.len() {
            let mut hashes = [0u64; Self::B];
            #[allow(clippy::needless_range_loop)]
            for j in 0..Self::B {
                hashes[j] = hash;
                hash = hash.rotate_left(1)
                    ^ wyhash_byte(*ptr.add(i + j))
                    ^ wyhash_byte(*ptr.add(i + j - k)).rotate_left(k as u32);
            }

            f(hashes, Self::B, i);

            i += Self::B;
        }

        let len = s.len() + 1 - i;

        if len > 0 {
            let mut hashes = [0u64; Self::B];
            #[allow(clippy::needless_range_loop)]
            for j in 0..len {
                hashes[j] = hash;
                hash = hash.rotate_left(1)
                    ^ wyhash_byte(*ptr.add(i + j))
                    ^ wyhash_byte(*ptr.add(i + j - k)).rotate_left(k as u32);
            }

            f(hashes, len, i);
        }
    }
}

impl SeedSearcher for GeneralSearcher {
    fn search(&self, text: &[u8], mut candidate_fn: impl FnMut(SeedMatch)) {
        unsafe {
            Self::get_hashes(text, self.k, |curr_hashes, len, text_i| {
                self.table
                    .batch_lookup::<{ Self::B }>(curr_hashes, len, text_i, &mut candidate_fn);
            });
        }
    }
}

struct HashToPatternIdx {
    filter: Filter,
    table: SeedTable,
    pattern_idxs: Vec<(u32, u32)>, // (pattern_idx, pattern_i)
}

impl HashToPatternIdx {
    pub fn new(mut hash_pattern_idxs: Vec<(u64, usize, usize)>) -> Self {
        assert!(hash_pattern_idxs.len() <= u32::MAX as usize);
        hash_pattern_idxs.sort_unstable();
        let filter = Filter::new(hash_pattern_idxs.iter().map(|(h, _, _)| *h));
        let table = SeedTable::new(&hash_pattern_idxs);

        Self {
            filter,
            table,
            pattern_idxs: hash_pattern_idxs
                .into_iter()
                .map(|(_, idx, i)| (idx as u32, i as u32))
                .collect(),
        }
    }

    #[inline(always)]
    pub unsafe fn batch_lookup<const N: usize>(
        &self,
        hashes: [u64; N],
        len: usize,
        text_i: usize,
        candidate_fn: &mut impl FnMut(SeedMatch),
    ) {
        let mut contains = 0u64;

        #[allow(clippy::needless_range_loop)]
        for i in 0..N {
            if self.filter.contains(hashes[i]) {
                contains |= 1 << i;
            }
        }

        contains &= (1 << len) - 1;

        while contains > 0 {
            let i = contains.trailing_zeros() as usize;
            contains &= contains - 1;

            let hash = hashes[i];
            let Some((start, end)) = self.table.get(hash) else {
                continue;
            };

            for &(pattern_idx, pattern_i) in &self.pattern_idxs[start as usize..end as usize] {
                candidate_fn(SeedMatch {
                    text_i: text_i + i,
                    pattern_idx: pattern_idx as usize,
                    pattern_i: pattern_i as usize,
                });
            }
        }
    }
}

#[cfg(not(feature = "seed-baseline"))]
struct SeedTable(HashTable<(u64, (u32, u32))>);

#[cfg(not(feature = "seed-baseline"))]
impl SeedTable {
    fn new(entries: &[(u64, usize, usize)]) -> Self {
        let mut table: HashTable<(u64, (u32, u32))> = HashTable::with_capacity(entries.len());
        for (i, (hash, _, _)) in entries.iter().enumerate() {
            if let Some((_, interval)) = table.find_mut(*hash, |(key, _)| key == hash) {
                interval.1 += 1;
            } else {
                table.insert_unique(*hash, (*hash, (i as u32, i as u32 + 1)), |(key, _)| *key);
            }
        }
        Self(table)
    }

    #[inline(always)]
    fn get(&self, hash: u64) -> Option<(u32, u32)> {
        self.0
            .find(hash, |(key, _)| *key == hash)
            .map(|(_, interval)| *interval)
    }
}

#[cfg(feature = "seed-baseline")]
struct SeedTable(FxHashMap<u64, (u32, u32)>);

#[cfg(feature = "seed-baseline")]
impl SeedTable {
    fn new(entries: &[(u64, usize, usize)]) -> Self {
        let mut table = FxHashMap::default();
        for (i, (hash, _, _)) in entries.iter().enumerate() {
            table.entry(*hash).or_insert((i as u32, i as u32)).1 += 1;
        }
        Self(table)
    }

    #[inline(always)]
    fn get(&self, hash: u64) -> Option<(u32, u32)> {
        self.0.get(&hash).copied()
    }
}

struct Filter {
    hashes: Vec<u16>,
    len: usize,
}

impl Filter {
    pub fn new(full_hashes: impl Iterator<Item = u64>) -> Self {
        let mut full_hashes = full_hashes.collect::<Vec<_>>();
        full_hashes.sort_unstable();
        full_hashes.dedup();
        let len = (full_hashes.len() * 2).next_power_of_two();
        let mut hashes = vec![0u16; len + 32];

        for h in full_hashes {
            Self::insert(&mut hashes, len, h);
        }

        Self { hashes, len }
    }

    fn insert(hashes: &mut [u16], len: usize, h: u64) {
        let mut idx = (h as usize) & (len - 1);

        loop {
            if idx >= hashes.len() {
                return;
            }
            if hashes[idx] == 0 {
                break;
            }
            idx += 1;
        }

        hashes[idx] = (((h >> 49) << 1) | 1) as u16;
    }

    #[inline(always)]
    pub unsafe fn contains(&self, hash: u64) -> bool {
        let idx = (hash as usize) & (self.len - 1);
        let hash_hi = (((hash >> 49) << 1) | 1) as u16;
        let mut match_mask = 0u64;
        let mut zero_mask = 0u64;

        cfg_if! {
            if #[cfg(any(target_feature = "avx2", all(target_arch = "x86_64", feature = "release-simd")))] {
                let hashes = _mm256_loadu_si256(self.hashes.as_ptr().add(idx) as _);
                let match_hash = _mm256_cmpeq_epi16(_mm256_set1_epi16(hash_hi as _), hashes);
                let match_zero = _mm256_cmpeq_epi16(hashes, _mm256_setzero_si256());
                match_mask |= _mm256_movemask_epi8(match_hash) as u32 as u64;
                zero_mask |= _mm256_movemask_epi8(match_zero) as u32 as u64;
            } else {
                for i in 0..16 {
                    if hash_hi == *self.hashes.as_ptr().add(idx + i) {
                        match_mask |= 1 << i;
                    }
                    if *self.hashes.as_ptr().add(idx + i) == 0 {
                        zero_mask |= 1 << i;
                    }
                }
            }
        }

        let first_matches_mask = (zero_mask & zero_mask.wrapping_neg()).wrapping_sub(1);
        match_mask &= first_matches_mask;
        zero_mask == 0 || match_mask > 0
    }
}

#[derive(Copy, Clone, Debug)]
pub struct SeedMatch {
    pub pattern_idx: usize,
    pub pattern_i: usize,
    pub text_i: usize,
}

#[repr(align(64))]
struct Aligned<const L: usize>([u8; L]);

cfg_if! {
    if #[cfg(any(target_feature = "avx2", all(target_arch = "x86_64", feature = "release-simd")))] {
        #[cfg(target_arch = "x86")]
        use std::arch::x86::*;
        #[cfg(target_arch = "x86_64")]
        use std::arch::x86_64::*;

        #[inline(always)]
        unsafe fn load_mask_avx2(len: usize) -> __m256i {
            static MASKS: [u8; 64] = {
                let mut m = [0u8; 64];
                let mut i = 0;

                while i < 32 {
                    m[i] = 0xFFu8;
                    i += 1;
                }

                m
            };
            _mm256_loadu_si256(MASKS.as_ptr().add(32 - len) as _)
        }

        #[inline(always)]
        unsafe fn read_rest_avx2(ptr: *const u8, len: usize) -> __m256i {
            let addr = ptr as usize;
            let start_page = addr >> 12;
            let end_page = (addr + 31) >> 12;

            if start_page == end_page {
                _mm256_loadu_si256(ptr as _)
            } else {
                let mut a = Aligned::<32>([0u8; 32]);
                std::ptr::copy_nonoverlapping(ptr, a.0.as_mut_ptr(), len);
                _mm256_load_si256(a.0.as_ptr() as _)
            }
        }
    }
}

#[inline(always)]
fn wyhash_byte(b: u8) -> u64 {
    #[cfg(feature = "seed-baseline")]
    {
        wyhash_byte_const(b)
    }
    #[cfg(not(feature = "seed-baseline"))]
    WYHASH_BYTE_LOOKUP[b as usize]
}

#[cfg(not(feature = "seed-baseline"))]
const WYHASH_BYTE_LOOKUP: [u64; 256] = {
    let mut lookup = [0u64; 256];
    let mut byte = 0;
    while byte < lookup.len() {
        lookup[byte] = wyhash_byte_const(byte as u8);
        byte += 1;
    }
    lookup
};

const fn wyhash_byte_const(byte: u8) -> u64 {
    let a = 0xa076_1d64_78bd_642fu64;
    let b = (byte as u64) ^ 0xe703_7ed1_a0b4_28dbu64;
    let c = (a as u128).wrapping_mul(b as u128);
    (c ^ (c >> 64)) as u64
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    fn direct_hash(kmer: &[u8]) -> u64 {
        kmer.iter()
            .fold(0u64, |hash, &byte| hash.rotate_left(1) ^ wyhash_byte(byte))
    }

    #[test]
    fn rolling_hash_matches_direct_recomputation_for_all_byte_values() {
        let sequence = (0u8..=255).collect::<Vec<_>>();
        for k in [1, 2, 7, 31, 64] {
            let mut observed = Vec::new();
            unsafe {
                GeneralSearcher::get_hashes(&sequence, k, |hashes, len, _| {
                    observed.extend_from_slice(&hashes[..len]);
                });
            }
            let expected = sequence.windows(k).map(direct_hash).collect::<Vec<_>>();
            assert_eq!(observed, expected, "rolling hash mismatch at k={k}");
        }
    }

    #[test]
    fn general_searcher_emits_every_exact_seed_candidate() {
        let patterns = [b"ACGTNacgt".as_slice(), b"NacgtXYZ".as_slice()];
        let text = b"--ACGTNacgt--NacgtXYZ--";
        let k = 4;
        let searcher = GeneralSearcher::new(patterns.iter().copied().enumerate(), k);
        let mut observed = BTreeSet::new();
        searcher.search(text, |seed| {
            observed.insert((seed.pattern_idx, seed.pattern_i - k, seed.text_i - k));
        });

        let mut expected = BTreeSet::new();
        for (pattern_idx, pattern) in patterns.iter().enumerate() {
            for (pattern_i, pattern_seed) in pattern.windows(k).enumerate() {
                for (text_i, text_seed) in text.windows(k).enumerate() {
                    if pattern_seed == text_seed {
                        expected.insert((pattern_idx, pattern_i, text_i));
                    }
                }
            }
        }
        assert!(expected.is_subset(&observed));
    }

    #[cfg(any(
        target_feature = "avx2",
        all(target_arch = "x86_64", feature = "release-simd")
    ))]
    #[test]
    fn small_searcher_emits_every_pattern_sharing_a_seed() {
        let patterns = [b"AC".as_slice(), b"AC".as_slice()];
        let searcher = SmallSearcher::<2>::new(patterns.iter().copied().enumerate()).unwrap();
        let mut observed = BTreeSet::new();
        searcher.search(b"AC", |seed| {
            observed.insert(seed.pattern_idx);
        });
        assert_eq!(observed, BTreeSet::from([0, 1]));
    }

    #[cfg(any(
        target_feature = "avx2",
        all(target_arch = "x86_64", feature = "release-simd")
    ))]
    #[test]
    fn small_searcher_emits_a_candidate_stored_in_bit_seven() {
        let patterns = [
            b"CC".as_slice(),
            b"CG".as_slice(),
            b"CT".as_slice(),
            b"GC".as_slice(),
            b"GG".as_slice(),
            b"GT".as_slice(),
            b"TC".as_slice(),
            b"AA".as_slice(),
        ];
        let searcher = SmallSearcher::<2>::new(patterns.iter().copied().enumerate()).unwrap();
        let mut observed = BTreeSet::new();
        searcher.search(b"AA", |seed| {
            observed.insert(seed.pattern_idx);
        });
        assert!(observed.contains(&7));
    }
}
