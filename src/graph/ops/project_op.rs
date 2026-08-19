use crate::graph::*;
use crate::inline_string::InlineString;

/// Apply a terminal output projection without materializing a concatenation.
///
/// The source labels must belong to the same FASTQ lane. Since projection
/// discards the lane's other mappings, this operation is intended for the end
/// of a transform graph immediately before output.
pub struct ProjectOp {
    required_names: Vec<LabelOrAttr>,
    str_type: StrType,
    labels: Vec<InlineString>,
}

impl ProjectOp {
    const NAME: &'static str = "ProjectOp";

    pub fn new(labels: impl IntoIterator<Item = Label>) -> Self {
        let labels = labels.into_iter().collect::<Vec<_>>();
        assert!(!labels.is_empty(), "ProjectOp requires at least one label");
        let str_type = labels[0].str_type;
        assert!(
            labels.iter().all(|label| label.str_type == str_type),
            "ProjectOp source labels must belong to one FASTQ lane"
        );
        let required_names = labels.iter().cloned().map(LabelOrAttr::Label).collect();
        let labels = labels.into_iter().map(|label| label.label).collect();
        Self {
            required_names,
            str_type,
            labels,
        }
    }
}

impl<T: Trace> GraphNode<T> for ProjectOp {
    fn run_inner(&self, mut reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        for read in &mut reads {
            read.project_whole(self.str_type, &self.labels)
                .map_err(|source| Error::NameError {
                    source,
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
