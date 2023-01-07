/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 Steven Walton
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.  
 *
 *
 * Class functions for an easier h5 read and write method.
 * Created by: Steven Walton
 * Email: walton.stevenj@gmail.com
 * Version: 1.2
 *
 * If you wish to contribute to this endeavour please email me or fork my branch from
 * https://github.com/stevenwalton/H5Easy
 * 
 */
#ifndef H5RW_H
#define H5RW_H

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <typeinfo>

#include "H5Cpp.h"

using namespace H5;

class WriteH5
{
   private:
      std::string variable;
      std::string filename;
   public:
      // sets our filename and our variable name
      void setFileName ( std::string name ) {filename = name;};
      void setVarName  ( std::string name ) {variable = name;};
      // Functions to be overloaded
      template<typename T>
      void writeData(const std::vector<T>&);
      template<typename T>
      void writeData(const std::vector<std::vector<T> >&);
      template<typename T>
      void writeData(const T&);

      void createGroup(std::string);
};

class LoadH5
{
   private:
      std::string variable;
      std::string filename;
   public:
      // sets our filename and our variable name
      void setFileName(std::string name) {filename = name;};
      void setVarName(std::string name) {variable = name;};
      // Read functions
      int getDataint() const;
      float getDatafloat() const;
      double getDatadouble() const;
      std::vector<int> getDataVint() const;
      std::vector<float> getDataVfloat() const;
      std::vector<double> getDataVDouble() const;
      std::vector<std::vector<int> > getData2Dint() const;
      std::vector<std::vector<float> > getData2Dfloat() const;
      std::vector<std::vector<double> > getData2Ddouble() const;
      // Return the size of the data
      // Note that for multi-dim arrays that it gets the total size and not the size of a single row.
      int getSize() const;

      // We now make a proxy class so that we can overload the return type and use a single
      // function to get data whether int or float. This could be made more advanced by 
      // adding more data types (such as double). 
      class Proxy
      {
         private:
            LoadH5 const* myOwner;
         public:
            Proxy( const LoadH5* owner ) : myOwner( owner ) {}
            operator int() const
            {
                return myOwner->getDataint();
            }
            operator float() const
            {
                return myOwner->getDatafloat();
            }
            operator double() const
            {
                return myOwner->getDatadouble();
            }
            operator std::vector<int>() const
            {
               return myOwner->getDataVint();
            }
            operator std::vector<float>() const
            {
               return myOwner->getDataVfloat();
            }
            operator std::vector<double>() const
            {
               return myOwner->getDataVDouble();
            }
            operator std::vector<std::vector<int> >() const
            {
               return myOwner->getData2Dint();
            }
            operator std::vector<std::vector<float> >() const
            {
               return myOwner->getData2Dfloat();
            }
            operator std::vector<std::vector<double> >() const
            {
               return myOwner->getData2Ddouble();
            }
      };
      // Here we use the Proxy class to have a single getData function
      Proxy getData() const {return Proxy(this);}
};

#endif
