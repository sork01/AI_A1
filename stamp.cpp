// *neode.onsave* setgo gcc -lstdc++ -lm stamp.cpp -o stamp && cat hmm4.sample | ./stamp
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <cmath>

#define putmat(m)         std::cout << m.numRows() << " " << m.numCols() << " ";\
        \
        for (unsigned int iii = 0; iii < m.numRows(); ++iii)\
        {\
            for (unsigned int jjj = 0; jjj < m.numCols(); ++jjj)\
            {\
                std::cout << m[iii][jjj] << " " ;\
            }\
        }\
        \
        std::cout << std::endl;

class vector
{
public:
    vector(): vector(1) {}
    
    vector(unsigned int size)
    {
        values.resize(size);
    }
    
    vector(std::vector<float> vals)
    {
        values = vals;
    }
    
    virtual ~vector() {}
    
    float& operator[](unsigned int i)
    {
        return values[i];
    }
    
    vector operator+(vector other)
    {
        vector result(size());
        
        for (unsigned int i = 0; i < size(); ++i)
        {
            result[i] = values[i] + other[i];
        }
        
        return result;
    }
    
    vector operator-(vector other)
    {
        vector result(size());
        
        for (unsigned int i = 0; i < size(); ++i)
        {
            result[i] = values[i] - other[i];
        }
        
        return result;
    }
    
    float dot(vector other)
    {
        float result;
        
        if (size() != other.size())
        {
            throw std::invalid_argument("(vector::dot) vector dimension mismatch");
        }
        else
        {
            for (unsigned int i = 0; i < size(); ++i)
            {
                result += values[i] * other[i];
            }
        }
        
        return result;
    }
    
    unsigned int size()
    {
        return values.size();
    }
    
    std::string toString()
    {
        std::ostringstream result;
        
        result << "(";
        
        for (int i = 0; i < values.size(); ++i)
        {
            if (i + 1 < values.size())
            {
                result << values[i] << ", ";
            }
            else
            {
                result << values[i];
            }
        }
        
        result << ")";
        
        return result.str();
    }

private:
    std::vector<float> values;
};

class matrix
{
public:
    matrix(): matrix(3,3) {}
    
    matrix(unsigned int rows, unsigned int cols)
    {
        m_rows = rows;
        m_cols = cols;
        
        for (unsigned int i = 0; i < rows; ++i)
        {
            row_vectors.emplace_back(vector(cols));
        }
    }
    
    vector& operator[](unsigned int i)
    {
        if (i >= m_rows)
        {
            throw std::invalid_argument("(matrix::operator[]) out of bounds");
        }
        else
        {
            return row_vectors[i];
        }
    }
    
    matrix operator*(matrix other)
    {
        if (m_cols != other.m_rows)
        {
            throw std::invalid_argument("(matrix::operator*) invalid dimensions");
        }
        else
        {
            matrix result(m_rows, other.m_cols);
            
            for (unsigned int i = 0; i < m_rows; ++i)
            {
                for (unsigned int j = 0; j < other.m_cols; ++j)
                {
                    // c[i,j] = a[row:i] dot b[col:j]
                    
                    for (unsigned int n = 0; n < m_cols; ++n)
                    {
                        result[i][j] += row_vectors[i][n] * other.row_vectors[n][j];
                    }
                }
            }
            
            return result;
        }
    }
    
    matrix operator-(matrix other)
    {
        if (m_rows != other.m_rows || m_cols != other.m_cols)
        {
            throw std::invalid_argument("(matrix::operator-) invalid dimensions");
        }
        else
        {
            matrix result(m_rows, m_cols);
            
            for (unsigned int i = 0; i < m_rows; ++i)
            {
                for (unsigned int j = 0; j < m_cols; ++j)
                {
                    result[i][j] = row_vectors[i][j] - other.row_vectors[i][j];
                }
            }
            
            return result;
        }
    }
    
    unsigned int numRows()
    {
        return m_rows;
    }
    
    unsigned int numCols()
    {
        return m_cols;
    }
    
    float sum()
    {
        float result = 0;
        
        for (unsigned int i = 0; i < m_rows; ++i)
        {
            for (unsigned int j = 0; j < m_cols; ++j)
            {
                result += row_vectors[i][j];
            }
        }
        
        return result;
    }
    
    std::string toString()
    {
        std::ostringstream result;
        
        for (vector v : row_vectors)
        {
            result << v.toString() << std::endl;
        }
        
        return result.str();
    }
    
    virtual ~matrix() {}

private:
    unsigned int m_rows;
    unsigned int m_cols;
    std::vector<vector> row_vectors;
};

class HMM
{
public:
    HMM() : HMM(matrix(), matrix(), matrix()) {}
    
    HMM(matrix A, matrix B, matrix pi) : m_A(A), m_B(B), m_pi(pi) {}
    
    matrix getNextEPD()
    {
        // question 2: what is the result of pi * A?
        //
        // pi * A gives an 1 x N vector representing the probability distribution P(X_1 | λ). the element
        // (pi*A)[i] gives the probability P(X_1 = x_i | λ) which is the probability that we started in
        // any state X_0 and then transitioned (in one step) to the state x_i
        
        // question 3: what is the result of (pi * A) * B?
        //
        // (pi * A) * B gives an 1 x N vector representing the probability distribution P(O_1 | λ).
        // the element ((pi*A)*B)[i] gives the probability P(O_1 = o_i | X_1, λ) which is the probability
        // of starting in any state X_0, transitioning (in one step) to X_1 and observing the emission
        // o_i
        
        
        return (m_pi*m_A)*m_B;
    }
    
    float getOSP()
    {
        computeAlpha();
        
        float result = 0;
        
        for (unsigned int i = 0; i < m_A.numRows(); ++i)
        {
            result += m_alpha[i][m_obseq.size()-1];
        }
        
        return result;
    }
    
    std::vector<int> getViterbi()
    {
        // question 4: why can we substitute O[1:t] = o[1:t] with O_t = o_t when conditioning on
        // X_t = x_i?
        //
        // this is probably because of "the markov property", i.e "given the present, the future
        // is independent of the past."
        
        unsigned int T = m_obseq.size();
        
        std::vector<int> result(T);
        
        matrix delta(m_A.numRows(), T);
        matrix delta_idx(m_A.numRows(), T);
        
        // question 5: how many values are stored in delta and delta_idx?
        //
        // in terms of memory, in our code these both store N x T values where N is the number of
        // states and T is the length of the sequence of observations. however, delta_idx does not
        // actually store any predecessors in its first column, so it would be fair to say that
        // delta_idx stores N x (T-1) values. maybe we should adapt the code to this fact.
        // delta (non-idx) stores a value in each and every cell and so is truly N x T.
        
        for (unsigned int i = 0; i < m_A.numRows(); ++i)
        {
            delta[i][0] = m_pi[0][i]*m_B[i][m_obseq[0]];
        }
        
        for (unsigned int t = 1; t < T; ++t)
        {
            for (unsigned int i = 0; i < m_A.numRows(); ++i)
            {
                float probmax = 0;
                unsigned int argmax = 0;
                
                for (unsigned int i2 = 0; i2 < m_A.numRows(); ++i2)
                {
                    float i2prob = delta[i2][t-1] * m_A[i2][i] * m_B[i][m_obseq[t]];
                    
                    if (i2prob > probmax)
                    {
                        probmax = i2prob;
                        argmax = i2;
                    }
                }
                
                delta[i][t] = probmax;
                delta_idx[i][t] = argmax;
            }
        }
        
        float probmax = 0;
        unsigned int argmax = 0;
        
        for (unsigned int i = 0; i < m_A.numRows(); ++i)
        {
            if (delta[i][T-1] > probmax)
            {
                probmax = delta[i][T-1];
                argmax = i;
            }
        }
        
        result[T-1] = argmax;
        unsigned int back = delta_idx[argmax][T-1];
        
        for (int t = T-2; t >= 0; --t)
        {
            // std::cout << "t mofo: " << t << std::endl << std::endl;
            result[t] = back;
            back = delta_idx[back][t];
        }
        
        return result;
    }
    
    void setObservationSequence(std::vector<int> seq)
    {
        // store the given sequence
        m_obseq = seq;
        
        m_need_alpha = true;
        m_need_beta = true;
        m_need_gamma = true;
        m_need_digamma = true;
        
        
        /*
        // number of states
        unsigned int N = m_A.numRows();
        
        // number of observations / time steps
        unsigned int T = seq.size();
        
        // for performance we dont want to actually compute all of these everytime a new sequence
        // is given (e.g if user just wants to do viterbi we dont want to compute digamma), but we
        // atleast reset them here.
        
        // alpha, beta and gamma is N x T
        m_alpha = matrix(N, T);
        m_beta = matrix(N, T);
        m_gamma = matrix(N, T);
        
        // digamma is N x N x T
        m_digamma = std::vector<matrix>(N);
        
        for (unsigned int i = 0; i < N; ++i)
        {
            m_digamma[i] = matrix(N, T);
        }
        */
    }
    
    bool train()
    {
        // std::cout << "training" << std::endl;
        
        computeDigamma();
        computeGamma();
        
        unsigned int N = m_A.numRows();
        unsigned int K = m_B.numCols();
        unsigned int T = m_obseq.size();
        
        std::vector<float> gamma_rowsum(N);
        
        for (unsigned int i = 0; i < N; ++i)
        {
            for (unsigned int t = 0; t < T-1; ++t)
            {
                gamma_rowsum[i] += m_gamma[i][t]; // ∑ t=1..(T-1) of γ_t(i) (2.35, 2.36)
            }
        }
        
        matrix newA(N, N);
        matrix newB(N, K);
        matrix newPi(1, N);
        
        for (unsigned int i = 0; i < N; ++i)
        {
            for (unsigned int j = 0; j < N; ++j)
            {
                for (unsigned int t = 0; t < T-1; ++t)
                {
                    newA[i][j] += m_digamma[t][i][j]; // 2.35 numerator
                }
                
                newA[i][j] /= gamma_rowsum[i]; // 2.35 denominator
            }
        }
        
        for (unsigned int k = 0; k < K; ++k)
        {
            for (unsigned int j = 0; j < N; ++j)
            {
                for (unsigned int t = 0; t < T-1; ++t)
                {
                    if (m_obseq[t] == k)
                    {
                        newB[j][k] += m_gamma[j][t]; // (2.36 numerator)
                    }
                }
                
                newB[j][k] /= gamma_rowsum[j]; // (2.36 denominator)
            }
        }
        
        for (unsigned int i = 0; i < N; ++i)
        {
            newPi[0][i] = m_gamma[i][0]; // (2.37)
        }
        
        
        
        // if (delta < min_delta)
        // {
            // std::cout << "aint doin shit" << std::endl;
            // return false;
        // }
        // else
        // {
            // std::cout << "practice makes perfect" << std::endl;
            
            m_A = newA;
            m_B = newB;
            m_pi = newPi;
            
            return true;
        // }
    }
    
    void stamp(int maxIters)
    {
        // 1. Initialization:
        // Select initial values for the matrices A, B and π, where π is 1 × N , while A = {a ij }
        // is N × N and B = {b j (k)} is N × M , and all three matrices are row-stochastic. If
        // known, use reasonable approximations for the matrix values, otherwise let π i ≈ 1/N
        // and a ij ≈ 1/N and b j (k) ≈ 1/M . Be sure that each row sums to 1 and the elements
        // of each matrix are not uniform. Then initialize
        // maxIters = maximum number of re-estimation iterations
        
        bool quit = false;
        
        int N = m_A.numRows();
        int M = m_B.numCols();
        int T = m_obseq.size();
        
        int iters = 0;
        float oldLogProb = -std::numeric_limits<float>::infinity();
        float logProb = 0;
        
        while (!quit)
        {
            matrix alpha(T, N);
            std::vector<float> c(T);
            matrix beta(T, N);
            matrix gamma(T, N);
            std::vector<matrix> digamma(T);
            
            for (unsigned int i = 0; i < T; ++i)
            {
                digamma[i] = matrix(N, N);
            }
            
            // 2. The α-pass
            // compute α 0 (i)
            // c 0 = 0
            c[0] = 0.0;
            // for i = 0 to N − 1
            for (unsigned int i = 0; i < N; ++i)
            {
                // α 0 (i) = π i b i (O 0 )
                alpha[0][i] = m_pi[0][i] * m_B[i][m_obseq[0]];
                // c 0 = c 0 + α 0 (i)
                
                c[0] += alpha[0][i];
                // next i
            }
            
            // scale the α 0 (i)
            // c 0 = 1/c 0
            c[0] = 1.0/c[0];
            
            // for i = 0 to N − 1
            for (unsigned int i = 0; i < N; ++i)
            {
                // α 0 (i) = c 0 α 0 (i)
                alpha[0][i] *= c[0];
                // next i
            }
            
            // compute α t (i)
            // for t = 1 to T − 1
            for (unsigned int t = 1; t < T; ++t)
            {
                // c t = 0
                c[t] = 0;
                
                // for i = 0 to N − 1
                for (unsigned int i = 0; i < N; ++i)
                {
                    // α t (i) = 0
                    alpha[t][i] = 0;
                    
                    // for j = 0 to N − 1
                    for (unsigned int j = 0; j < N; ++j)
                    {
                        // α t (i) = α t (i) + α t−1 (j)a ji
                        alpha[t][i] += alpha[t-1][j] * m_A[j][i];
                        // next j
                    }
                    
                    // α t (i) = α t (i)b i (O t )
                    alpha[t][i] *= m_B[i][m_obseq[t]];
                    
                    // c t = c t + α t (i)
                    c[t] += alpha[t][i];
                    // next i
                }
                
                // scale α t (i)
                // c t = 1/c t
                c[t] = 1.0/c[t];
                
                // for i = 0 to N − 1
                for (unsigned int i = 0; i < N; ++i)
                {
                    // α t (i) = c t α t (i)
                    alpha[t][i] *= c[t];
                    // next i
                }
                
                // next t
            }
            
            // 3. The β-pass
            // Let β T −1 (i) = 1, scaled by c T −1
            // for i = 0 to N − 1
            for (unsigned int i = 0; i < N; ++i)
            {
                // β T −1 (i) = c T −1
                beta[T-1][i] = c[T-1];
                // next i
            }
            
            // β-pass
            // for t = T − 2 to 0 by − 1
            for (int t = T-2; t >= 0; --t)
            {
                // for i = 0 to N − 1
                for (unsigned int i = 0; i < N; ++i)
                {
                    // β t (i) = 0
                    beta[t][i] = 0;
                    
                    // for j = 0 to N − 1
                    for (unsigned int j = 0; j < N; ++j)
                    {
                        // β t (i) = β t (i) + a ij b j (O t+1 )β t+1 (j)
                        beta[t][i] += m_A[i][j] * m_B[j][m_obseq[t+1]] * beta[t+1][j];
                        // next j
                    }
                    
                    // scale β t (i) with same scale factor as α t (i)
                    // β t (i) = c t β t (i)
                    beta[t][i] *= c[t];
                    // next i
                }
                // next t
            }
            
            // 4. Compute γ t (i, j) and γ t (i)
            // for t = 0 to T − 2
            for (unsigned int t = 0; t < T-1; ++t)
            {
                // denom = 0
                float denom = 0.0;
                
                // for i = 0 to N − 1
                for (unsigned int i = 0; i < N; ++i)
                {
                    // for j = 0 to N − 1
                    for (unsigned int j = 0; j < N; ++j)
                    {
                        // denom = denom + α t (i)a ij b j (O t+1 )β t+1 (j)
                        denom += alpha[t][i] * m_A[i][j] * m_B[j][m_obseq[t+1]] * beta[t+1][j];
                        
                        // next j
                    }
                    
                    // next i
                }
                
                // for i = 0 to N − 1
                for (unsigned int i = 0; i < N; ++i)
                {
                    // γ t (i) = 0
                    gamma[t][i] = 0;
                    
                    // for j = 0 to N − 1
                    for (unsigned int j = 0; j < N; ++j)
                    {
                        // γ t (i, j) = (α t (i)a ij b j (O t+1 )β t+1 (j))/denom
                        digamma[t][i][j] = (alpha[t][i] * m_A[i][j] * m_B[j][m_obseq[t+1]] * beta[t+1][j]) / denom;
                        // γ t (i) = γ t (i) + γ t (i, j)
                        gamma[t][i] += digamma[t][i][j];
                        // next j
                    }
                    
                    // next i
                }
                
                // next t
            }
            
            // Special case for γ T −1 (i)
            // denom = 0
            float denom = 0;
            
            // for i = 0 to N − 1
            for (unsigned int i = 0; i < N; ++i)
            {
                // denom = denom + α T −1 (i)
                denom += alpha[T-1][i];
                // next i
            }
            
            // for i = 0 to N − 1
            for (unsigned int i = 0; i < N; ++i)
            {
                // γ T −1 (i) = α T −1 (i)/denom
                gamma[T-1][i] = alpha[T-1][i] / denom;
                // next i
            }
            
            // 5. Re-estimate A, B and π
            // re-estimate π
            // for i = 0 to N − 1
            for (unsigned int i = 0; i < N; ++i)
            {
            // π i = γ 0 (i)
            // next i
            }
            
            // re-estimate A
            // for i = 0 to N − 1
            for (unsigned int i = 0; i < N; ++i)
            {
                // for j = 0 to N − 1
                for (unsigned int j = 0; j < N; ++j)
                {
                    // numer = 0
                    float numer = 0;
                    // denom = 0
                    float denom = 0;
                    
                    // for t = 0 to T − 2
                    for (unsigned int t = 0; t < T-1; ++t)
                    {
                        // numer = numer + γ t (i, j)
                        numer += digamma[t][i][j];
                        
                        // denom = denom + γ t (i)
                        denom += gamma[t][i];
                        // next t
                    }
                    
                    // a ij = numer/denom
                    m_A[i][j] = numer/denom;
                    
                    // next j
                }
                
                // next i
            }
            
            // re-estimate B
            // for i = 0 to N − 1
            for (unsigned int i = 0; i < N; ++i)
            {
                // for j = 0 to M − 1
                for (unsigned int j = 0; j < M; ++j)
                {
                    // numer = 0
                    float numer = 0;
                    
                    // denom = 0
                    float denom = 0;
                    
                    // for t = 0 to T − 1
                    for (unsigned int t = 0; t < T; ++t)
                    {
                        // if(O t == j) then
                        if (m_obseq[t] == j)
                        {
                            // numer = numer + γ t (i)
                            numer += gamma[t][i];
                            
                            // end if
                        }
                        
                        // denom = denom + γ t (i)
                        denom += gamma[t][i];
                        
                        // next t
                    }
                    
                    // b i (j) = numer/denom
                    m_B[i][j] = numer/denom;
                    
                    // next j
                }
                
                // next i
            }
            
            // 6. Compute log[P (O | λ)]
            // logProb = 0
            logProb = 0;
            
            // for i = 0 to T − 1
            for (unsigned int i = 0; i < T; ++i)
            {
                // logProb = logProb + log(c i )
                logProb += log10(c[i]);
                
                // next i
            }
            
            // logProb = −logProb
            logProb = -logProb;
            
            // 7. To iterate or not to iterate, that is the question. . .
            // iters = iters + 1
            iters += 1;
            
            // if (iters < maxIters and logProb > oldLogProb) then
            if (iters < maxIters && logProb > oldLogProb)
            {
                // oldLogProb = logProb
                
                // printf("hey im all better: %.50f vs %.50f\n", logProb, oldLogProb);
                
                oldLogProb = logProb;
                
                // goto 2
            }// else
            else
            {
                // printf("well that sucked: %.50f vs %.50f\n", logProb, oldLogProb);
                quit = true;
                // output λ = (π, A, B)
                // end if
            }
        }
    }
    
    matrix getTransitionMatrix()
    {
        return m_A;
    }
    
    matrix getEmissionMatrix()
    {
        return m_B;
    }

private:
    matrix m_A;
    matrix m_B;
    matrix m_pi;
    
    bool m_need_alpha = false;
    bool m_need_beta = false;
    bool m_need_gamma = false;
    bool m_need_digamma = false;
    
    matrix m_alpha;
    matrix m_alpha_scaled;
    matrix m_beta;
    matrix m_gamma;
    
    std::vector<matrix> m_digamma;
    
    std::vector<int> m_obseq;
    
    std::vector<float> c;
    
    void computeAlpha()
    {
        printf("computing alpha\n");
        /*
        if (!m_need_alpha)
        {
            // dont recompute alpha if not necessary
            return;
        }
        
        m_need_alpha = false;
        */
        
        unsigned int N = m_A.numRows();
        unsigned int T = m_obseq.size();
        
        printf("%d %d\n", N, T);
        
        // alpha is N x T
        m_alpha = matrix(N, T);
        
        for (unsigned int i = 0; i < N; ++i)
        {
            m_alpha[i][0] = m_pi[0][i] * m_B[i][m_obseq[0]];
        }
        
        for (unsigned int t = 1; t < T; ++t)
        {
            for (unsigned int i = 0; i < N; ++i)
            {
                m_alpha[i][t] = 0;
                
                for (unsigned int i2 = 0; i2 < N; ++i2)
                {
                    m_alpha[i][t] += m_alpha[i2][t-1] * m_A[i2][i] * m_B[i][m_obseq[t]];
                }
            }
        }
    }
    
    void computeBeta()
    {
        // we need c[]
        computeAlpha();
        
        if (!m_need_beta)
        {
            // dont recompute beta if not necessary
            return;
        }
        
        m_need_beta = false;
        
        unsigned int N = m_A.numRows();
        unsigned int T = m_obseq.size();
        
        // beta is N x T
        m_beta = matrix(N, T);
        
        // final column of beta is given and constant
        for (unsigned int i = 0; i < N; ++i)
        {
            // m_beta[i][T-1] = 1;
            m_beta[i][T-1] = c[T-1]; // i.e 1 * c[T-1]
        }
        
        for (int t = T-2; t >= 0; --t) // step backwards in time
        {
            for (unsigned int i = 0; i < N; ++i) // row to set in beta
            {
                for (unsigned int i2 = 0; i2 < N; ++i2) // sum rows of next column in beta
                {
                    m_beta[i][t] += m_beta[i2][t+1] * m_B[i2][m_obseq[t+1]] * m_A[i][i2];
                }
                
                m_beta[i][t] * c[t];
            }
        }
    }
    
    void computeDigamma()
    {
        if (!m_need_digamma)
        {
            // dont recompute digamma if not necessary
            return;
        }
        
        m_need_digamma = false;
        
        computeAlpha();
        computeBeta();
        
        unsigned int T = m_obseq.size();
        unsigned int N = m_A.numRows();
        
        // digamma is T x N x N (T-1 ?????)
        m_digamma = std::vector<matrix>(T);
        
        for (unsigned int t = 0; t < T; ++t)
        {
            m_digamma[t] = matrix(N, N);
        }
        
        // now we have a set of T empty N x N matrices
        
        // the denominator looks constant . sum(k=1..N, alpha_T(k))
        
        // so compute it first
        
        float denom = 0;
        
        for (unsigned int k = 0; k < N; ++k)
        {
            // std::cout << "lol m_alpha[" << k << "][" << (T-1) << "] is " << m_alpha[k][T-1] << std::endl;
            denom += m_alpha_scaled[k][T-1];
        }
        
        // std::cout << "lol denom is " << denom << std::endl;
        // exit(0);
        
        // now a stack of loops
        
        for (unsigned int t = 0; t < T-1; ++t)
        {
            for (unsigned int i = 0; i < N; ++i)
            {
                for (unsigned int j = 0; j < N; ++j)
                {
                    m_digamma[t][i][j] = m_alpha_scaled[i][t] * m_A[i][j] * m_B[j][m_obseq[t+1]] * m_beta[j][t+1];
                    m_digamma[t][i][j] /= denom;
                    
                    // question 6: why do we need the division here?
                    //
                    // the denominator is the sum of the final column of the alpha matrix
                    // this is the same sum as is returned by getOSP, i.e the total probability
                    // for a given sequence of observations independent of what states and what order
                    // of states generated the sequence.
                    //
                    // in a probability-theoretic sense, dividing by this number suggests that we are
                    // constraining ourselves to a world in which the sequence of observations has
                    // taken place. this can be seen symbolically in equation 2.31 where we are
                    // looking at the conditional probability given that the sequence of observations
                    // is fixed (O_[1:T] = o_[1:T]), which corresponds to figuring out a probability
                    // in the sample space constrained to the world where O_[1:T] = o_[1:T] holds,
                    // which is precisely what is accomplished by the division. i think.
                }
            }
        }
    }
    
    void computeGamma()
    {
        if (!m_need_gamma)
        {
            // dont recompute gamma if not necessary
            return;
        }
        
        m_need_gamma = false;
        
        computeDigamma();
        
        // for (unsigned int tired = 0; tired < m_obseq.size(); ++tired)
        // {
            // putmat(m_digamma[tired]);
        // }
        
        // std::cout << "bajs" << std::endl;
        // exit(0);
        
        unsigned int N = m_A.numRows();
        unsigned int T = m_obseq.size();
        
        // gamma is N x T (T-1? TO SAVE LOTS OF SPACE)
        m_gamma = matrix(N, T);
        
        for (unsigned int t = 0; t < T-1; ++t)
        {
            for (unsigned int i = 0; i < N; ++i)
            {
                for (unsigned int j = 0; j < N; ++j)
                {
                    m_gamma[i][t] += m_digamma[t][i][j];
                }
            }
        }
    }
};

class labhelper
{
public:
    static matrix parseMatrix()
    {
        unsigned int rows;
        unsigned int cols;
        
        std::cin >> rows >> cols;
        
        matrix result(rows, cols);
        
        for (unsigned int i = 0; i < rows; ++i)
        {
            for (unsigned int j = 0; j < cols; ++j)
            {
                std::cin >> result[i][j];
            }
        }
        
        return result;
    }
    
    static std::vector<int> parseIntSequence()
    {
        unsigned int n;
        std::cin >> n;
        std::vector<int> result(n);
        
        for (unsigned int i = 0; i < n; ++i)
        {
            std::cin >> result[i];
        }
        
        return result;
    }
    
    static void printMatrix(matrix m)
    {
        std::cout << m.numRows() << " " << m.numCols() << " ";
        
        for (unsigned int i = 0; i < m.numRows(); ++i)
        {
            for (unsigned int j = 0; j < m.numCols(); ++j)
            {
                // std::cout << m[i][j] << " " ;
                printf("%.6f ", m[i][j]);
            }
        }
        
        std::cout << std::endl;
    }
    
    static void printIntSequence(std::vector<int> seq)
    {
        for (unsigned int i = 0; i < seq.size(); ++i)
        {
            std::cout << seq[i] << " ";
        }
        
        std::cout << std::endl;
    }
};

int _main()
{
    matrix A = labhelper::parseMatrix();
    matrix B = labhelper::parseMatrix();
    matrix pi = labhelper::parseMatrix();
    std::vector<int> o = labhelper::parseIntSequence();
    
    HMM model(A,B,pi);
    
    // std::cout << model.getOSP(o) << std::endl;
    model.setObservationSequence(o);
    std::cout << model.getOSP() << std::endl;    
}

int main()
{
    matrix A = labhelper::parseMatrix();
    matrix B = labhelper::parseMatrix();
    matrix pi = labhelper::parseMatrix();
    std::vector<int> o = labhelper::parseIntSequence();
    
    HMM model(A,B,pi);
    
    // labhelper::printMatrix(model.getTransitionMatrix());
    // labhelper::printMatrix(model.getEmissionMatrix());
    // exit(0);
    
    model.setObservationSequence(o);
    
    // for (int i = 0; i < 1000; ++i)
    // {
        // std::cout << "[" << i << "] sequence prob: " << model.getOSP() << std::endl;
        // model.train();
    // }
    
    model.stamp(100);
    labhelper::printMatrix(model.getTransitionMatrix());
    labhelper::printMatrix(model.getEmissionMatrix());
    
    // matrix m1(2,2);
    // m1[0][0] = 1;
    // m1[0][1] = 2;
    // m1[1][0] = 3;
    // m1[1][1] = 4;
    
    // matrix m2 = m1;
    // labhelper::printMatrix(m2);
    
    // model.train(0.1);
    
    // std::cout << model.getOSP() << std::endl;
    // labhelper::printIntSequence(model.getViterbi());
}
