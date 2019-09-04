/*
MIT License

Copyright (c) 2019 Tucker Haydon

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <cstdlib>

// RAII buffer for storing a 2D floating point array
class SmartBuffer2D
{
public:
    SmartBuffer2D(
        const size_t rows,
        const size_t cols)
        : rows_(rows),
          cols_(cols)
    {
        this->data_ = new double *[this->rows_];
        for (size_t row = 0; row < this->rows_; ++row)
        {
            this->data_[row] = new double[this->cols_];
        }
    }

    ~SmartBuffer2D()
    {
        for (size_t row = 0; row < this->rows_; ++row)
        {
            delete[] this->data_[row];
        }
        delete[] this->data_;
    }

    double **Get() const
    {
        return this->data_;
    }

    size_t Rows() const
    {
        return this->rows_;
    }

    size_t Cols() const
    {
        return this->cols_;
    }

private:
    double **data_;
    size_t rows_;
    size_t cols_;
};

// RAII buffer for storing a 1D floating point array
class SmartBuffer1D
{
public:
    SmartBuffer1D(const size_t rows)
        : rows_(rows)
    {
        this->data_ = new double[this->rows_];
    }

    ~SmartBuffer1D()
    {
        delete[] this->data_;
    }

    double *Get() const
    {
        return this->data_;
    }

    size_t Rows() const
    {
        return this->rows_;
    }

private:
    double *data_;
    size_t rows_;
};