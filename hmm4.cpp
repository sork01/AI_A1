// *neode.onsave* setgo gcc -g -lstdc++ -lm hmm4.cpp -o hmm4 && cat hmm4.sample | ./hmm4
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
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
    
    /**
     * Computes the probability distribution of emissions for the next step of
     * the current HMM.
     * 
     * EPD = Emission Probability Distribution.
     *
     * @return An 1xN matrix representing the vector P(O[t+1] | λ)
     */
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
    
    /**
     * Computes the overall probability of the current HMM ever generating the given
     * observation sequence regardless of which states are involved.
     * 
     * OSP = Observation Sequence Probability.
     * 
     * @return The probability in decimal
     */
    float _getOSP()
    {
        computeAlpha();
        
        // return pow(10, getLogOSP());
        
        unsigned int N = m_A.numRows();
        unsigned int T = m_obseq.size();
        
        float result = 0;
        
        for (unsigned int i = 0; i < N; ++i)
        {
            result += m_alpha[i][T-1];
        }
        
        // return log10(result);
        return result;
    }
    
    float getOSP() // log(0.090276) ~~ -1.04443
    {
        return pow(10, getLogOSP());
    }
    
    float getLogOSP()
    {
        if (m_obseq.size() < 1)
        {
            throw std::logic_error("(HMM::getLogOSP) empty observation sequence");
        }
        
        computeAlpha();
        
        float result = 0.0f;
        
        for (unsigned int i = 0; i < m_obseq.size(); ++i)
        {
            // printf("c[%d] = %.6f\n", i, c[i]);
            result += log10(c[i]);
        }
        
        return -result;
    }
    
    /**
     * getViterbi computes the most likely series of states that the current HMM would
     * go through in order to generate the given sequence of observations. The computation
     * uses the Viterbi algorithm.
     */
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
    
    /**
     *
     */
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
        
        float oldlogprob = getLogOSP();
        printf("oldlogprob = %.6f\n", oldlogprob);
        
        
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
        
        putmat(newPi)
        exit(0);
        
        /*
        if (delta < min_delta)
        {
            std::cout << "aint doin shit" << std::endl;
            return false;
        }
        else
        {
            std::cout << "practice makes perfect" << std::endl;
            
            m_A = newA;
            m_B = newB;
            m_pi = newPi;
            
            return true;
        }
        */
        
        HMM newModel = HMM(newA, newB, newPi);
        newModel.setObservationSequence(m_obseq);
        
        float logprob = newModel.getLogOSP();
        
        
        printf("newlogprob = %.6f\n", logprob);
        
        if (abs(logprob - oldlogprob) > 1.0)
        {
            printf("practice makes perfect\n");
            
            m_A = newA;
            m_B = newB;
            m_pi = newPi;
            
            return true;
        }
        else
        {
            printf("aint doin shit\n");
            return false;
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
    
    matrix getInitialStatePD()
    {
        return m_pi;
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
    matrix m_beta;
    matrix m_gamma;
    std::vector<matrix> m_digamma;
    
    std::vector<int> m_obseq;
    std::vector<float> c;
    
    void computeAlpha()
    {
        if (!m_need_alpha)
        {
            // dont recompute alpha if not necessary
            return;
        }
        
        m_need_alpha = false;
        
        unsigned int N = m_A.numRows();
        unsigned int T = m_obseq.size();
        
        c = std::vector<float>(T);
        
        // alpha is N x T
        m_alpha = matrix(N, T);
        
        for (unsigned int i = 0; i < T; ++i)
        {
            c[i] = 0;
        }
        
        for (unsigned int i = 0; i < N; ++i)
        {
            m_alpha[i][0] = m_pi[0][i] * m_B[i][m_obseq[0]];
            c[0] += m_alpha[i][0];
        }
        
        c[0] = 1.0/c[0];
        
        // printf("c[%d] = %.6f\n", 0, c[0]);
        
        for (unsigned int i = 0; i < N; ++i)
        {
            m_alpha[i][0] *= c[0];
        }
        
        for (unsigned int t = 1; t < T; ++t)
        {
            for (unsigned int i = 0; i < N; ++i)
            {
                for (unsigned int i2 = 0; i2 < N; ++i2)
                {
                    m_alpha[i][t] += m_alpha[i2][t-1] * m_A[i2][i] * m_B[i][m_obseq[t]];
                    // c[t] += m_alpha[i][t];
                }
                
                c[t] += m_alpha[i][t];
                
                // printf("alpha[t=%d,i=%d] = %.6f\n", t, i, m_alpha[i][t]);
                // printf("%.6f ", t, i, m_alpha[i][t]);
            }
            
            // printf("\n");
            
            c[t] = 1.0/c[t];
            
            // printf("c[%d] = %.6f\n", t, c[t]);
            
            for (unsigned int i = 0; i < N; ++i)
            {
                m_alpha[i][t] *= c[t];
                // printf("%.6f ", t, i, m_alpha[i][t]);
            }
            
            // printf("\n");
        }
        
        // computeAlphaScaled();
    }
    
    // void computeAlphaScaled()
    // {
        // unsigned int N = m_A.numRows();
        // unsigned int T = m_obseq.size();
        
        // m_alpha_scaled = m_alpha;
        
        // c = std::vector<float>(T);
        
        // for (unsigned int i = 0; i < N; ++i)
        // {
            // c[0] += m_alpha[i][0];
        // }
        
        // c[0] = 1.0/c[0];
        
        // for (unsigned int i = 0; i < N; ++i)
        // {
            // m_alpha_scaled[i][0] = m_alpha[i][0] * c[0];
        // }
        
        // for (unsigned int t = 1; t < T; ++t)
        // {
            // for (unsigned int i = 0; i < N; ++i)
            // {
                // for (unsigned int i2 = 0; i2 < N; ++i2)
                // {
                    // m_alpha_scaled[i][t] += m_alpha_scaled[i2][t-1] * m_A[i2][i] * m_B[i][m_obseq[t]];
                // }
                
                // c[t] += m_alpha_scaled[i][t];
            // }
            
            // c[t] = 1.0/c[t];
            
            // for (unsigned int i = 0; i < N; ++i)
            // {
                // m_alpha_scaled[i][t] *= c[t];
            // }
        // }
    // }
    
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
        
        for (unsigned int t = 0; t < T-1; ++t)
        {
            float denom = 0.0;
            
            for (unsigned int i = 0; i < N; ++i)
            {
                for (unsigned int j = 0; j < N; ++j)
                {
                    denom += m_alpha[i][t] * m_A[i][j] * m_B[j][m_obseq[t+1]] * m_beta[j][t+1];
                }
            }
            
            for (unsigned int i = 0; i < N; ++i)
            {
                for (unsigned int j = 0; j < N; ++j)
                {
                    m_digamma[t][i][j] = m_alpha[i][t] * m_A[i][j] * m_B[j][m_obseq[t+1]] * m_beta[j][t+1];
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
                std::cout << m[i][j] << " " ;
            }
        }
        
        std::cout << std::endl;
    }
    
    static void printMatrix2(matrix m)
    {
        printf("%d %d ", m.numRows(), m.numCols());
        
        for (unsigned int i = 0; i < m.numRows(); ++i)
        {
            for (unsigned int j = 0; j < m.numCols(); ++j)
            {
                printf("%.6f ", m[i][j]);
            }
        }
        
        printf("\n");
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

// hmm1
int hmm1main()
{
    matrix A = labhelper::parseMatrix();
    matrix B = labhelper::parseMatrix();
    matrix pi = labhelper::parseMatrix();
    std::vector<int> o = labhelper::parseIntSequence();
    
    HMM model(A,B,pi);
    
    model.setObservationSequence(o);
    printf("%.6f\n", model.getOSP());
}

// hmm2
int hmm2main()
{
    
    matrix A = labhelper::parseMatrix();
    matrix B = labhelper::parseMatrix();
    matrix pi = labhelper::parseMatrix();
    std::vector<int> o = labhelper::parseIntSequence();
    
    HMM model(A,B,pi);
    
    model.setObservationSequence(o);
    labhelper::printIntSequence(model.getViterbi());
}

// hmm4
int main()
{
    matrix A = labhelper::parseMatrix();
    matrix B = labhelper::parseMatrix();
    matrix pi = labhelper::parseMatrix();
    std::vector<int> o = labhelper::parseIntSequence();
    
    HMM model(A,B,pi);
    model.setObservationSequence(o);
    
    for (int i = 0; i < 1; ++i)
    {
        model.train();
    }
    
    labhelper::printMatrix2(model.getTransitionMatrix());
}
