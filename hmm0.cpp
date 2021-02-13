// *neode.onsave* setgo gcc -lstdc++ hmm0.cpp -o hmm0 && cat hmm0.sample | ./hmm0
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

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
    
    unsigned int numRows()
    {
        return m_rows;
    }
    
    unsigned int numCols()
    {
        return m_cols;
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
        return (m_pi*m_A)*m_B;
    }

private:
    matrix m_A;
    matrix m_B;
    matrix m_pi;
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
};

int main()
{
    
    matrix A = labhelper::parseMatrix();
    matrix B = labhelper::parseMatrix();
    matrix pi = labhelper::parseMatrix();
    
    HMM model(A,B,pi);
    
    labhelper::printMatrix(model.getNextEPD());
}
