# First discrepancy: fixed-base are not faster than variable-base MSMs (when using variable-time)
cargo bench -- "Multiscalar multiplications/Variable-time variable-base multiscalar multiplication/1024$"
cargo bench -- "Multiscalar multiplications/Variable-time fixed-base multiscalar multiplication/1024$"

# Second discrepancy: constant-time fixed-base scalar mul is faster than variable-time fixed-base scalar mul
# Possible explanation: the variable-time fixed-base scalar mul is implemented as an MSM; some overhead from there?
cargo bench -- "Constant-time fixed-base scalar mul"
cargo bench -- "Multiscalar multiplications/Variable-time fixed-base multiscalar multiplication/1$"
