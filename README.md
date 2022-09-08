# OrderEfficientMatChainMul
Order Efficient Matrix Chain Multiplication, implemented in C

The order in which one executes a series of matrix multiplications can drastically impact the number of operations required to find the result. For example, if one wants to multiply three matrices, A, B, and C, then multiplying BC and only then A(BC) can save an (arbitrarily) large number of operations and execute much more quickly compared to first multiplying (AB) and then (AB)C. Thus, it is faster first to calculate the order which involves the fewest multiplications and then execute in that order.
Wikipedia goes more in-depth on this https://en.wikipedia.org/wiki/Matrix_chain_multiplication

For those who want faster matrix multiplication, enjoy!
