use regex::bytes::*;

use thread_local::*;

use crate::graph::*;
use crate::inline_string::*;

pub struct MatchRegexOp {
    required_names: Vec<LabelOrAttr>,
    produced_names: Vec<LabelOrAttr>,
    label: Label,
    attr: Option<Attr>,
    regex: Regex,
    capture_names: Vec<InlineString>,
    regex_local: ThreadLocal<Regex>,
}

impl MatchRegexOp {
    const NAME: &'static str = "MatchRegexOp";

    /// Match a regex pattern in an interval.
    ///
    /// If named capture groups are used, then intervals are automatically created at the match
    /// locations, labeled by the names specified in the regex.
    ///
    /// The transform expression must have one input label and one output attribute.
    ///
    /// Example `transform_expr`: `tr!(seq1.* -> seq1.*.matched)`.
    /// This will match the regex pattern two `seq1.*` and set `seq1.*.matched` to a boolean
    /// indicating whether the regex matches.
    pub fn new(transform_expr: TransformExpr, regex: &str) -> Self {
        Self::try_new(transform_expr, regex).unwrap_or_else(|error| panic!("{error}"))
    }

    pub fn try_new(transform_expr: TransformExpr, regex: &str) -> Result<Self> {
        transform_expr.try_check_size(1, 1, Self::NAME)?;
        transform_expr.try_check_same_str_type(Self::NAME)?;

        let label = transform_expr.before(0);
        let attr = transform_expr.try_after_attr(0, Self::NAME)?;
        let regex = Regex::new(regex).map_err(|error| Error::InvalidOperation {
            operation: Self::NAME,
            reason: format!("invalid regular expression: {error}"),
        })?;
        let capture_names = regex
            .capture_names()
            .filter_map(|name| name.map(|name| InlineString::new(name.as_bytes())))
            .collect::<Vec<_>>();
        let mut produced_names = attr
            .iter()
            .cloned()
            .map(LabelOrAttr::Attr)
            .collect::<Vec<_>>();
        produced_names.extend(capture_names.iter().map(|&name| {
            LabelOrAttr::Label(Label {
                str_type: label.str_type,
                label: name,
            })
        }));

        Ok(Self {
            required_names: vec![label.clone().into()],
            produced_names,
            label,
            attr,
            regex,
            capture_names,
            regex_local: ThreadLocal::new(),
        })
    }
}

impl<T: Trace> GraphNode<T> for MatchRegexOp {
    fn produced_names(&self) -> Option<&[LabelOrAttr]> {
        Some(&self.produced_names)
    }

    fn effects_are_complete(&self) -> bool {
        true
    }

    fn mutation_kind(&self) -> MutationKind {
        MutationKind::Metadata
    }

    fn rejection_behavior(&self) -> RejectionBehavior {
        RejectionBehavior::Never
    }

    fn cost_class(&self) -> CostClass {
        CostClass::Search
    }

    fn run_inner(&self, mut reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        let regex = self.regex_local.get_or(|| self.regex.clone());

        for read in &mut reads {
            let mut new_mappings = Vec::new();

            let string = read
                .substring(self.label.str_type, self.label.label)
                .map_err(|e| Error::NameError {
                    source: e,
                    read: read.clone(),
                    context: Self::NAME,
                })?;
            let matched;

            match regex.captures(string) {
                Some(caps) => {
                    matched = true;

                    new_mappings.extend(self.capture_names.iter().filter_map(|&name| {
                        caps.name(name.as_str()).map(|m| (name, m.start(), m.len()))
                    }));
                }
                None => matched = false,
            }

            let str_mappings = read.str_mappings_mut(self.label.str_type).unwrap();
            let offset = str_mappings.mapping(self.label.label).unwrap().start;

            for (label, start, len) in new_mappings.drain(..) {
                str_mappings.add_mapping(Some(label), offset + start, len);
            }

            if let Some(attr) = &self.attr {
                let target = match read.data_mut(attr.str_type, attr.label, attr.attr) {
                    Ok(target) => target,
                    Err(source) => {
                        return Err(Error::NameError {
                            source,
                            read: read.clone(),
                            context: Self::NAME,
                        });
                    }
                };
                *target = Data::Bool(matched);
            }
        }

        Ok((Some(reads), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
