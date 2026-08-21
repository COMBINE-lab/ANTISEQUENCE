use crate::graph::*;

pub struct NullOutputOp;

impl Default for NullOutputOp {
    fn default() -> Self {
        Self::new()
    }
}

impl NullOutputOp {
    const NAME: &'static str = "NullOutputOp";

    pub fn new() -> Self {
        Self
    }
}

impl<T: Trace> GraphNode<T> for NullOutputOp {
    fn stage(&self) -> NodeStage {
        NodeStage::Output
    }

    fn produced_names(&self) -> Option<&[LabelOrAttr]> {
        Some(&[])
    }

    fn mutation_kind(&self) -> MutationKind {
        MutationKind::None
    }

    fn rejection_behavior(&self) -> RejectionBehavior {
        RejectionBehavior::Never
    }

    fn cost_class(&self) -> CostClass {
        CostClass::Constant
    }

    fn supports_prepared_output(&self) -> bool {
        true
    }

    fn has_explicit_name_observation(&self) -> bool {
        true
    }

    fn supports_direct_projection(&self) -> bool {
        true
    }

    fn produces_prepared_output(&self) -> bool {
        false
    }

    fn prepare_output(
        &self,
        _reads: &[Read],
        _recycled: Option<PreparedOutput>,
    ) -> Result<PreparedOutput> {
        Ok(PreparedOutput::Passthrough)
    }

    fn prepare_output_projected(
        &self,
        _reads: &[Read],
        _projections: &[DirectReadProjection<'_>],
        _recycled: Option<PreparedOutput>,
    ) -> Result<PreparedOutput> {
        Ok(PreparedOutput::Passthrough)
    }

    fn commit_output(&self, prepared: &mut PreparedOutput) -> Result<()> {
        if matches!(prepared, PreparedOutput::Passthrough) {
            Ok(())
        } else {
            Err(Error::InvalidPipelineGraph(
                "NullOutputOp received an incompatible prepared payload".to_owned(),
            ))
        }
    }

    fn run_inner(&self, reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        Ok((Some(reads), false))
    }

    fn required_names(&self) -> &[LabelOrAttr] {
        &[]
    }

    fn name(&self) -> &'static str {
        Self::NAME
    }
}
