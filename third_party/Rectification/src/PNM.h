//-*-c++-*-

#ifndef _PNM_H
#define _PNM_H

#include "Rectify/Ops.h"

#include <iostream>
#include <string>

class PNMLoadSave {

private:

  // Skip over spaces and comments; temp is the current file character
  static void SkipSpaces(std::ifstream& file, signed char& temp)
  {
    while (temp == ' ' || temp == '\t' || temp == '\n' || temp == '#')
    {
      if (file.eof())
        return;
      if (temp == '#') // skip this line:
        while ((temp != '\n') && (temp != EOF))
          temp=file.get();
      // skip this `whitespace' byte:
      temp=file.get();
    }
  }

  // Get an integer from the file stream; temp is the current file character
  static int ReadInteger(std::ifstream& file, signed char& temp)
  {
    int n = 0;
    bool neg=false;
    while (((temp >= '0') && (temp <= '9')) || (temp == '-'))
    {
      if (file.eof())
        return 0;
      if (temp=='-')
        neg=true;
      n *= 10; n += (temp - '0');
      temp=file.get();
    }
    if (neg)
      n=-n;
    return n;
  }

  struct pgm_image_header
  {
    char magic_number[2]; // magic number
    int height; // number of scan lines
    int width; // scan line length
    int maxval;
  };

public:

  /// Load the image in the given filename. If succeeds, image is filled and returns true. Returns false otherwise
  template <class DataT>
    static bool LoadImage(const std::string &filename, Image::Image<DataT> &image) {
      // assume filename already verified
      std::ifstream file(filename.c_str(), std::ios::in | std::ios::binary);
      if (!file) {
        std::cerr << "Couldn't open file '" << filename << "' for reading image" << std::endl;
        exit(1);
      }

      unsigned int xsize, ysize, bitspp, magicnumber;
      bool colour;

      colour=false; // PGM is greyscale by default
      signed char temp;

      /** READ THE HEADER **/

      //This method assumes the PGM header is in the following format
      //P[25]
      //#comments
      //#up to any amount of comments (including no comments)
      //width height
      //maxval
      file.seekg(0, std::ios::beg);
      //Because PGM Image files can contain an unlimited number of
      //comments, this method of reading in the header is infeasible
      //It was left in to read in the magic number
      pgm_image_header header;
      file.read((char*)(&header), sizeof(header));

      if ((header.magic_number[0]!='P') || (!file)) {
        std::cerr << "Couldn't open image file! " << std::endl;
        exit(1);
      }
      magicnumber=header.magic_number[1]-'0';
      if (magicnumber!=2 && magicnumber!=3 && magicnumber!=5 && magicnumber!=6) {
        std::cerr << "Couldn't open pnm image file because found bad magic number: only 2,3,5,6 recognised" << std::endl;
        exit(1);
      }

      file.seekg(3, std::ios::beg);
      temp=file.get();

      //Skip over spaces and comments
      SkipSpaces(file,temp);
      //Read in Width
      xsize=(unsigned int)ReadInteger(file,temp);
      //Skip over spaces and comments
      SkipSpaces(file,temp);
      //Read in Height
      ysize=ReadInteger(file,temp);
      //Skip over spaces and comments
      SkipSpaces(file,temp);
      //Read in Maxval
      unsigned int maxval=ReadInteger(file,temp);

      if (temp == ' ' || temp == '\t' || temp == '\n')
        temp=file.get();

      // Step the file back one.
      file.seekg(-1,std::ios::cur);

      if (maxval==0)
        bitspp=8;
      else
      {
        bitspp=0;
        while (maxval>0)
        {
          bitspp++; maxval=maxval>>1;
        }
      }

      /** DONE READING THE HEADER **/

      if (magicnumber==3 || magicnumber==6) {
        // Colour so must have 3 of those samples per pixel
        bitspp=bitspp*3;
        colour=true;
      }

      image.resize(xsize, ysize, colour);

      // Every image should have 1D access unless the image class is written by a retard
      // However, I'm going to assume retard level and simply assume a 1D block for the image at (0,0)
      // so don't add unnecesary functions to Image.h and freak those retards out
      if (magicnumber==2)
      {
        DataT TempStore;
        DataT *oiter=&(image(0,0)), *containerend=oiter+(xsize*ysize);

        while (!file.eof() && oiter!=containerend) {
          // Problems with my gcc unless copy to a temp - probably optimisations or something
          file >> TempStore;
          *oiter++=TempStore;
        }
      }
      if (magicnumber==3)
      {
        DataT R,G,B;

        DataT *oiter=&(image(0,0)), *containerend=oiter+(xsize*ysize);
        while (!file.eof() && oiter!=containerend) {
          file >> R >> G >> B;
          *oiter++=(R<<16)+(G<<8)+B;
        }
      }
      if (magicnumber==5)
      {
        DataT *oiter=&(image(0,0)), *containerend=oiter+(xsize*ysize);
        while (!file.eof() && oiter!=containerend)
          *oiter++=file.get();
      }
      if (magicnumber==6)
      {
        DataT *oiter=&(image(0,0)), *containerend=oiter+(xsize*ysize);
        while (!file.eof() && oiter!=containerend)
          *oiter++=(file.get()<<16)+(file.get()<<8)+file.get();
      }

      if (!(file.good()) && !(file.eof()))
        return false;

      return true;
    }

    // Save an image to the given file.
    template <class DataT>
    static void SaveImage(const std::string &filename, const Image::Image<DataT> &image) {
      std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary);
      if (!file) {
        std::cerr << "Couldn't open file '" << filename << "' for writing image" << std::endl;
        exit(1);
      }

      if (!image.isColour()) {
        DataT maxval=0;
        // Get bits per pixel for container based on max value if greyscale
        for (size_t x=0;x<image.xsize();++x) {
          for (size_t y=0;y<image.ysize();++y) {
            if (image(x,y)>maxval)
              maxval=image(x,y);
          }
        }
        // Get bits per pixel based on maximum value in the image
        unsigned int bpp=0;
        while (maxval>0) {
          ++bpp;
          maxval=maxval>>1;
        }

        // Write magic number.
        if (bpp==8)
          file << "P5" << std::endl;
        else
          file << "P2" << std::endl;

        // Write comments
        file << "# Created by some nice code" << std::endl;

        // Write Height and Width
        file << (unsigned int)image.xsize() << " " << (unsigned int)image.ysize() << std::endl;

        // Write max val
        if (bpp<32)
          file << (1<<bpp)-1 << std::endl;
        else
          file << UINT_MAX << std::endl;

        if (bpp==8) {
          char *buffer=new char[image.xsize()*image.ysize()];
          char *iter=buffer;
          for (size_t y=0;y<image.ysize();++y)
            for (size_t x=0;x<image.xsize();++x)
              *iter++=image(x,y);
          file.write(buffer, ((std::streamsize)(image.xsize()*image.ysize())));
          delete[] buffer;
        }
        else
        {
          for (size_t y=0;y<image.ysize();++y)
            for (size_t x=0;x<image.xsize();++x)
              file << image(x,y) << std::endl;
        }
      }
      else
      {

        file << "P6" << std::endl;

        // Write comments
        file << "# Created by some nice code" << std::endl;

        // Write Height and Width
        file << (unsigned int)image.xsize() << " " << (unsigned int)image.ysize() << std::endl;

        // Write max val
        file << "255" << std::endl;

        for (size_t y=0;y<image.ysize();++y)
          for (size_t x=0;x<image.xsize();++x) {
            DataT imval=image(x,y);
            file.put((imval & (255<<16))>>16);
            file.put((imval & (255<<8))>>8);
            file.put(imval & 255);
          }
      }

      if (!file) {
        std::cerr << "Filesystem error whilst writing image: " << filename << std::endl;
        exit(1);
      }

      return;

    }

};

#endif
