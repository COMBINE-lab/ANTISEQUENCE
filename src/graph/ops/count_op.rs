use parking_lot::Mutex;
use thread_local::ThreadLocal;

use crate::graph::*;

pub struct CountOp {
    required_names: Vec<LabelOrAttr>,
    selector_exprs: Vec<Expr>,
    counts: ThreadLocal<Mutex<Vec<usize>>>,
}

impl CountOp {
    const NAME: &'static str = "CountOp";

    /// For each selector expression, count the number of reads where the expression evaluates to
    /// true.
    pub fn new(selector_exprs: impl IntoIterator<Item = impl Into<Expr>>) -> Self {
        let selector_exprs = selector_exprs
            .into_iter()
            .map(|e| e.into())
            .collect::<Vec<_>>();
        let required_names = selector_exprs
            .iter()
            .flat_map(|n| n.required_names())
            .collect();
        Self {
            required_names,
            selector_exprs,
            counts: ThreadLocal::new(),
        }
    }

    /// Returns the counts.
    pub fn counts(&self) -> Vec<usize> {
        let mut totals = vec![0usize; self.selector_exprs.len()];
        for local in self.counts.iter() {
            for (total, count) in totals.iter_mut().zip(local.lock().iter()) {
                *total += count;
            }
        }
        totals
    }
}

impl<T: Trace> GraphNode<T> for CountOp {
    fn produced_names(&self) -> Option<&[LabelOrAttr]> {
        Some(&[])
    }

    fn mutation_kind(&self) -> MutationKind {
        MutationKind::None
    }

    fn rejection_behavior(&self) -> RejectionBehavior {
        RejectionBehavior::Never
    }

    fn run_inner(&self, reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        let counts = self
            .counts
            .get_or(|| Mutex::new(vec![0usize; self.selector_exprs.len()]));
        let mut counts = counts.lock();
        for read in &reads {
            for (count, selector) in counts.iter_mut().zip(&self.selector_exprs) {
                if selector.eval_bool(read).map_err(|e| Error::NameError {
                    source: e,
                    read: read.clone(),
                    context: Self::NAME,
                })? {
                    *count += 1;
                }
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
