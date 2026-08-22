use antisequence::graph::{GraphBuilder, InputFastqOp, MissingInputPolicy, NullOutputOp};
use antisequence::trace::NoTrace;
use std::io::Cursor;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut builder =
        GraphBuilder::<NoTrace>::new().with_missing_input_policy(MissingInputPolicy::Error);
    builder.add(InputFastqOp::from_reader(Cursor::new(
        b"@read\nACGT\n+\nIIII\n".as_slice(),
    ))?);
    builder.add(NullOutputOp::new());
    builder.compile()?.run()?;
    Ok(())
}
