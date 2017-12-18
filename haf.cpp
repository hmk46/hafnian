#include <iterator>
#include <cstdint>


long long hafnian(const long long mat[], const uint32_t rank)
{
	uint32_t N = rank/2;
	long long B[3][rank][rank][rank+1] = {0};
	long long h[rank+1] = {};
	h[N] = 1;
	for(int i = rank; i--; ) {
		for(int j = rank; j--; ) {
			B[0][i][j][0] = mat[rank*i+j];
		}
	}
	int32_t p, q; // loop variables
	int32_t z, w; // polynomial multiplication loop variables
	long long co; // accumulator for polynomial multiplication
	long long *B0p, *B1p; // convenience pointers for computation
	uint8_t i_r, i_w; // read and write indices that get swapped after each iteration
	uint32_t deg, off; // degree of current iteration polynomials and indexing offset
	long long *beta; // convenience pointer to polynomial accumulator
	uint32_t dim; // dimension of sub-matrix
	int8_t parity = 1-2*(N&1);
	int8_t sign; // sign of subset polynomial
	uint64_t s=1<<N, ss, b; // loop variables for powerset iteration
	while(ss = --s) {
		sign = parity;
		dim = rank-2;
		i_r = 0;
		i_w = 1;
		long long g[rank+1] = {};
		g[0] = 1;
		deg = 0;
		off = 0;
		// iterate through least significant bits (members) of the subset
		while(ss) {
			// adjust the offset and dimension based on the count of trailing off bits
			for(b = 1; !(b&ss); b<<=1) dim-=2, off+=2;
			sign *= -1;
			beta = B[i_r][off][1+off];
			// multiply g by 1+x beta
			// loop through coefficients backwards to multiply in-place
			for(z = deg; z>=0; z--) {
				co = g[z];
				for(w = 1; w <= deg+1 && w <= rank-z; w++) {
					g[z+w] += co * beta[w-1];
				}
			}
			for(p = 0; p < dim; p++) {
				B0p = B[i_r][off][p+2+off];
				B1p = B[i_r][1+off][p+2+off];
				for(q = p+1; q < dim; q++) {
					B[i_w][p][q][0] = B[i_r][p+2+off][q+2+off][0];
					// combined truncated polynomial multiplication and addition
					for(z = rank-1; z>=0; z--) {
						B[i_w][p][q][z+1] = B[i_r][p+2+off][q+2+off][z+1];
						for(w = 0; w <= deg && w <= rank-1-z; w++) {
							B[i_w][p][q][w+z+1] += B0p[w]*B[i_r][1+off][q+2+off][z]+B[i_r][off][q+2+off][w]*B1p[z];
						}
					}
				}
			}
			deg = std::min(rank, 2*deg+1);
			// cycle read and write indices
			i_r = i_w;
			i_w ^= 3;
			// reset read offset and shrink submatrix
			off = 0;
			dim -= 2;
			// divide out the bit vector by 1+ the power of the LSB (effectively clearing it)
			ss /= (ss^(ss&(ss-1)))<<1;
		}
		// add x^N g to h
		for(z = N; z <= rank; z++) h[z] += sign*g[z-N];
	}
	return h[rank];
}

int main(int argc, char **argv) {
	uint64_t rank = atoi(argv[1]);
	long long mat[rank*rank];
	for(uint64_t i = 0; i < rank; i++) {
		for(uint64_t j = 0; j < rank; j++) {
			mat[i*rank+j] = (i==j ? 0.0 : 1.0);
		}
	}
	long long h = hafnian(mat, rank);
	printf("%llu\n", h);
	return 0;
}
