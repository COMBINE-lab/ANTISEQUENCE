use crate::graph::*;
use crate::read::ReadProjectionPart;

/// One component of a terminal output projection.
#[derive(Clone, Debug)]
pub enum ProjectPart {
    Label(Label),
    Literal(Vec<u8>),
}

impl From<Label> for ProjectPart {
    fn from(label: Label) -> Self {
        Self::Label(label)
    }
}

impl ProjectPart {
    pub fn literal(bytes: impl Into<Vec<u8>>) -> Self {
        Self::Literal(bytes.into())
    }
}

/// Apply a terminal output projection without materializing a concatenation.
///
/// The source labels must belong to the same FASTQ lane. Since projection
/// discards the lane's other mappings, this operation is intended for the end
/// of a transform graph immediately before output.
pub struct ProjectOp {
    required_names: Vec<LabelOrAttr>,
    str_type: StrType,
    parts: Vec<ProjectPart>,
}

impl ProjectOp {
    const NAME: &'static str = "ProjectOp";

    pub fn new(labels: impl IntoIterator<Item = Label>) -> Self {
        let labels = labels.into_iter().collect::<Vec<_>>();
        assert!(!labels.is_empty(), "ProjectOp requires at least one label");
        let str_type = labels[0].str_type;
        Self::with_parts(str_type, labels.into_iter().map(ProjectPart::Label))
    }

    /// Create a projection that may interleave source labels and fixed bytes.
    pub fn with_parts(str_type: StrType, parts: impl IntoIterator<Item = ProjectPart>) -> Self {
        let parts = parts.into_iter().collect::<Vec<_>>();
        assert!(!parts.is_empty(), "ProjectOp requires at least one part");
        assert!(
            parts.iter().all(
                |part| !matches!(part, ProjectPart::Label(label) if label.str_type != str_type)
            ),
            "ProjectOp source labels must belong to one FASTQ lane"
        );
        let required_names = parts
            .iter()
            .filter_map(|part| match part {
                ProjectPart::Label(label) => Some(LabelOrAttr::Label(label.clone())),
                ProjectPart::Literal(_) => None,
            })
            .collect();
        Self {
            required_names,
            str_type,
            parts,
        }
    }
}

impl<T: Trace> GraphNode<T> for ProjectOp {
    fn run_inner(&self, mut reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        let parts = self
            .parts
            .iter()
            .map(|part| match part {
                ProjectPart::Label(label) => ReadProjectionPart::Mapping(label.label),
                ProjectPart::Literal(bytes) => ReadProjectionPart::Literal(bytes),
            })
            .collect::<Vec<_>>();
        for read in &mut reads {
            read.project_whole_with_literals(self.str_type, &parts)
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
