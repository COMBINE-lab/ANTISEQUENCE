use rand::distributions::Bernoulli;
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;

use crate::graph::*;

pub struct BernoulliOp {
    required_names: Vec<LabelOrAttr>,
    produced_names: Vec<LabelOrAttr>,
    attr: Attr,
    bernoulli: Bernoulli,
    seed: u64,
}

impl BernoulliOp {
    const NAME: &'static str = "BernoulliOp";

    /// Set the attribute `attr` to a sampled boolean from a Bernoulli distribution
    /// with probability `prob` of true.
    ///
    /// This is fully deterministic for a chosen seed and ordering of reads, even with
    /// multiple threads.
    pub fn new(attr: Attr, prob: f64, seed: u32) -> Self {
        Self::try_new(attr, prob, seed).unwrap_or_else(|error| panic!("{error}"))
    }

    pub fn try_new(attr: Attr, prob: f64, seed: u32) -> Result<Self> {
        let required_names = vec![LabelOrAttr::Label(Label {
            str_type: attr.str_type,
            label: attr.label,
        })];
        let produced_names = vec![LabelOrAttr::Attr(attr.clone())];
        let bernoulli = Bernoulli::new(prob).map_err(|error| Error::InvalidOperation {
            operation: Self::NAME,
            reason: error.to_string(),
        })?;
        Ok(Self {
            required_names,
            produced_names,
            attr,
            bernoulli,
            seed: seed as u64,
        })
    }
}

impl<T: Trace> GraphNode<T> for BernoulliOp {
    fn produced_names(&self) -> Option<&[LabelOrAttr]> {
        Some(&self.produced_names)
    }

    fn mutation_kind(&self) -> MutationKind {
        MutationKind::Metadata
    }

    fn rejection_behavior(&self) -> RejectionBehavior {
        RejectionBehavior::Never
    }

    fn cost_class(&self) -> CostClass {
        CostClass::Constant
    }

    fn run_inner(&self, mut reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        for read in &mut reads {
            // use the index of the read in the seed for determinism when multithreading
            let seed = (self.seed << 32).wrapping_add(read.first_idx() as u64);
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
            let rand_bool = self.bernoulli.sample(&mut rng);

            let target = match read.data_mut(self.attr.str_type, self.attr.label, self.attr.attr) {
                Ok(target) => target,
                Err(source) => {
                    return Err(Error::NameError {
                        source,
                        read: read.clone(),
                        context: Self::NAME,
                    });
                }
            };
            *target = Data::Bool(rand_bool);
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
