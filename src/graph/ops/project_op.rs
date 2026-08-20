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
        Self::try_new(labels).unwrap_or_else(|error| panic!("{error}"))
    }

    pub fn try_new(labels: impl IntoIterator<Item = Label>) -> Result<Self> {
        let labels = labels.into_iter().collect::<Vec<_>>();
        let Some(first) = labels.first() else {
            return Err(Error::InvalidOperation {
                operation: Self::NAME,
                reason: "at least one source label is required".to_owned(),
            });
        };
        let str_type = first.str_type;
        Self::try_with_parts(str_type, labels.into_iter().map(ProjectPart::Label))
    }

    /// Create a projection that may interleave source labels and fixed bytes.
    pub fn with_parts(str_type: StrType, parts: impl IntoIterator<Item = ProjectPart>) -> Self {
        Self::try_with_parts(str_type, parts).unwrap_or_else(|error| panic!("{error}"))
    }

    /// Fallible projection constructor with literal support.
    pub fn try_with_parts(
        str_type: StrType,
        parts: impl IntoIterator<Item = ProjectPart>,
    ) -> Result<Self> {
        let parts = parts.into_iter().collect::<Vec<_>>();
        if parts.is_empty() {
            return Err(Error::InvalidOperation {
                operation: Self::NAME,
                reason: "at least one projection part is required".to_owned(),
            });
        }
        if parts
            .iter()
            .any(|part| matches!(part, ProjectPart::Label(label) if label.str_type != str_type))
        {
            return Err(Error::InvalidOperation {
                operation: Self::NAME,
                reason: "all source labels must belong to the projected FASTQ lane".to_owned(),
            });
        }
        let required_names = parts
            .iter()
            .filter_map(|part| match part {
                ProjectPart::Label(label) => Some(LabelOrAttr::Label(label.clone())),
                ProjectPart::Literal(_) => None,
            })
            .collect();
        Ok(Self {
            required_names,
            str_type,
            parts,
        })
    }
}

impl<T: Trace> GraphNode<T> for ProjectOp {
    fn direct_read_projection(&self) -> Option<DirectReadProjection<'_>> {
        Some(DirectReadProjection {
            str_type: self.str_type,
            parts: &self.parts,
        })
    }

    fn produced_names(&self) -> Option<&[LabelOrAttr]> {
        Some(&[])
    }

    fn mutation_kind(&self) -> MutationKind {
        MutationKind::Record
    }

    fn rejection_behavior(&self) -> RejectionBehavior {
        RejectionBehavior::Never
    }

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
