use crate::graph::*;

pub struct RetainOp {
    required_names: Vec<LabelOrAttr>,
    selector_expr: Expr,
    constant_result: Option<bool>,
}

impl RetainOp {
    const NAME: &'static str = "RetainOp";

    /// Retain only the reads where the selector expression evaluates to true and discard the rest.
    pub fn new(selector_expr: impl Into<Expr>) -> Self {
        let mut selector_expr = selector_expr.into();
        let constant_result = selector_expr
            .optimize()
            .then(|| selector_expr.eval_bool(&Read::new()).ok())
            .flatten();
        Self {
            required_names: selector_expr.required_names(),
            selector_expr,
            constant_result,
        }
    }
}

impl<T: Trace> GraphNode<T> for RetainOp {
    fn produced_names(&self) -> Option<&[LabelOrAttr]> {
        Some(&[])
    }

    fn effects_are_complete(&self) -> bool {
        true
    }

    fn mutation_kind(&self) -> MutationKind {
        MutationKind::None
    }

    fn is_semantic_noop(&self) -> bool {
        self.constant_result == Some(true)
    }

    fn is_selective_filter(&self) -> bool {
        self.constant_result != Some(true)
    }

    fn is_infallible_selective_filter(&self) -> bool {
        self.constant_result == Some(false)
    }

    fn run_inner(&self, mut reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        let mut error = None;
        reads.retain(|read| {
            match self.selector_expr.eval_bool(read) {
                Ok(keep) => keep,
                Err(e) => {
                    error = Some(Error::NameError {
                        source: e,
                        read: read.clone(),
                        context: Self::NAME,
                    });
                    false // drop read if error? or stop?
                }
            }
        });

        if let Some(e) = error {
            return Err(e);
        }

        if reads.is_empty() {
            Ok((None, false))
        } else {
            Ok((Some(reads), false))
        }
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &self.required_names
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
