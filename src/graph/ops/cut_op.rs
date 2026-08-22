use crate::graph::*;

pub struct CutOp {
    required_names: Vec<LabelOrAttr>,
    produced_names: Vec<LabelOrAttr>,
    cut_label: Label,
    new_label1: Option<Label>,
    new_label2: Option<Label>,
    cut_idx: Expr,
}

impl CutOp {
    const NAME: &'static str = "CutOp";

    /// Cut a labeled interval at the specified index to create two new intervals.
    ///
    /// The transform expression must have one input label and two output labels.
    ///
    /// Example `transform_expr`: `tr!(seq1.* -> seq1.left, seq1.right)`.
    pub fn new(transform_expr: TransformExpr, cut_idx: impl Into<Expr>) -> Self {
        let cut_idx = cut_idx.into();
        transform_expr.check_size(1, 2, Self::NAME);
        transform_expr.check_same_str_type(Self::NAME);
        let mut required_names = cut_idx.required_names();
        required_names.push(transform_expr.before(0).into());
        let new_label1 = transform_expr.after_label(0, Self::NAME);
        let new_label2 = transform_expr.after_label(1, Self::NAME);
        let produced_names = new_label1
            .iter()
            .chain(new_label2.iter())
            .cloned()
            .map(LabelOrAttr::Label)
            .collect();

        Self {
            required_names,
            produced_names,
            cut_label: transform_expr.before(0),
            new_label1,
            new_label2,
            cut_idx,
        }
    }
}

impl<T: Trace> GraphNode<T> for CutOp {
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

    fn run_inner(&self, mut reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        for read in &mut reads {
            let cut_idx = self.cut_idx.eval_int(read).map_err(|e| Error::NameError {
                source: e,
                read: read.clone(),
                context: Self::NAME,
            })?;

            read.cut(
                self.cut_label.str_type,
                self.cut_label.label,
                self.new_label1.as_ref().map(|l| l.label),
                self.new_label2.as_ref().map(|l| l.label),
                cut_idx,
            )
            .map_err(|e| Error::NameError {
                source: e,
                read: read.clone(),
                context: Self::NAME,
            })?;
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
