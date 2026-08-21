use crate::graph::*;

pub struct SetOp {
    required_names: Vec<LabelOrAttr>,
    label_or_attr: LabelOrAttr,
    expr: Expr,
}

impl SetOp {
    const NAME: &'static str = "SetOp";

    /// Set a labeled interval or attribute to the result of an expression.
    ///
    /// The expression must return a byte string if a labeled interval is being set.
    ///
    /// To generate the quality scores when setting intervals that have corresponding quality
    /// scores, references to intervals in the expression are directly substituted with the
    /// corresponding quality scores of the intervals. For references to byte strings without
    /// quality scores, a sequence of `I`s is used as the quality scores in the expression.
    /// *This naive substitution may lead to unexpected results for complex expressions!*
    ///
    /// If a label is set, then its interval and all other intersecting intervals will be adjusted accordingly
    /// for any shortening or lengthening.
    pub fn new(label_or_attr: impl Into<LabelOrAttr>, expr: impl Into<Expr>) -> Self {
        let label_or_attr = label_or_attr.into();
        let mut expr = expr.into();
        expr.optimize();
        let mut required_names = expr.required_names();
        match &label_or_attr {
            LabelOrAttr::Label(_) => required_names.push(label_or_attr.clone()),
            LabelOrAttr::Attr(a) => required_names.push(LabelOrAttr::Label(Label {
                str_type: a.str_type,
                label: a.label,
            })),
            LabelOrAttr::RecordAttr(_) | LabelOrAttr::LaneAttr(_) => {}
        }

        Self {
            required_names,
            label_or_attr,
            expr,
        }
    }
}

impl<T: Trace> GraphNode<T> for SetOp {
    fn produced_names(&self) -> Option<&[LabelOrAttr]> {
        Some(std::slice::from_ref(&self.label_or_attr))
    }

    fn mutation_kind(&self) -> MutationKind {
        match &self.label_or_attr {
            LabelOrAttr::Label(_) => MutationKind::Sequence,
            LabelOrAttr::Attr(_) | LabelOrAttr::RecordAttr(_) | LabelOrAttr::LaneAttr(_) => {
                MutationKind::Metadata
            }
        }
    }

    fn rejection_behavior(&self) -> RejectionBehavior {
        RejectionBehavior::Never
    }

    fn run_inner(&self, mut reads: Vec<Read>) -> Result<(Option<Vec<Read>>, bool)> {
        for read in &mut reads {
            match &self.label_or_attr {
                LabelOrAttr::Label(label) => {
                    let new_bytes = self
                        .expr
                        .eval_bytes(read, false)
                        .map_err(|e| Error::NameError {
                            source: e,
                            read: read.clone(),
                            context: Self::NAME,
                        })?
                        .into_owned();

                    let str_mappings =
                        read.str_mappings(label.str_type)
                            .ok_or_else(|| Error::NameError {
                                source: NameError::NotInRead(Name::StrType(label.str_type)),
                                read: read.clone(),
                                context: Self::NAME,
                            })?;

                    if str_mappings.qual().is_some() {
                        let new_qual = self
                            .expr
                            .eval_bytes(read, true)
                            .map_err(|e| Error::NameError {
                                source: e,
                                read: read.clone(),
                                context: Self::NAME,
                            })?
                            .into_owned();

                        read.set(label.str_type, label.label, &new_bytes, Some(&new_qual))
                            .map_err(|e| Error::NameError {
                                source: e,
                                read: read.clone(),
                                context: Self::NAME,
                            })?;
                    } else {
                        read.set(label.str_type, label.label, &new_bytes, None)
                            .map_err(|e| Error::NameError {
                                source: e,
                                read: read.clone(),
                                context: Self::NAME,
                            })?;
                    }
                }
                LabelOrAttr::Attr(attr) => {
                    let new_val = self.expr.eval(read, false).map_err(|e| Error::NameError {
                        source: e,
                        read: read.clone(),
                        context: Self::NAME,
                    })?;
                    let new_val: Data = new_val.into();

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
                    *target = new_val;
                }
                LabelOrAttr::RecordAttr(attr) => {
                    let new_val: Data = self
                        .expr
                        .eval(read, false)
                        .map_err(|source| Error::NameError {
                            source,
                            read: read.clone(),
                            context: Self::NAME,
                        })?
                        .into();
                    *read.record_data_mut(attr.attr) = new_val;
                }
                LabelOrAttr::LaneAttr(attr) => {
                    let new_val: Data = self
                        .expr
                        .eval(read, false)
                        .map_err(|source| Error::NameError {
                            source,
                            read: read.clone(),
                            context: Self::NAME,
                        })?
                        .into();
                    *read.lane_data_mut(attr.lane, attr.attr) = new_val;
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
