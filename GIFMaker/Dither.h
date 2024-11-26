#pragma once
/*
MIT License

Copyright (c) 2024 Chris Lomont

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

// Color dither
// Chris Lomont, 2023

#include <algorithm> // max
#include <cmath>     // ceil, log2
#include <cstdint>   // uint8_t

class Dither
{

    // Riemersma dither
    // based on Hilbert curve from
    // https://www.compuphase.com/riemer.htm
    // error propagation completely redone, seems better
    // want video stable dither https://surma.dev/things/ditherpunk/


    // floating point color for accuracy
    struct RGB
    {
        RGB(const uint8_t* rgb) : RGB(rgb[0], rgb[1], rgb[2]) {}
        RGB(double r = 0, double g = 0, double b = 0) : r(r), g(g), b(b) {}
        double r{ 0 }, g{ 0 }, b{ 0 };
        [[nodiscard]] double lenSq() const { return r * r + g * g + b * b; }
        void Draw(uint8_t* rgb) const
        {
            rgb[0] = static_cast<uint8_t>(r);
            rgb[1] = static_cast<uint8_t>(g);
            rgb[2] = static_cast<uint8_t>(b);
        }
        friend RGB operator+(const RGB& lhs, const RGB& rhs) { return RGB(lhs.r + rhs.r, lhs.g + rhs.g, lhs.b + rhs.b); }
        friend RGB operator-(const RGB& lhs, const RGB& rhs) { return RGB(lhs.r - rhs.r, lhs.g - rhs.g, lhs.b - rhs.b); }
        friend RGB operator*(double lhs, const RGB& rhs) { return RGB(lhs * rhs.r, lhs * rhs.g, lhs * rhs.b); }
    };
    
    int xpos{}, ypos{};
    int width{}, height{};
    uint8_t* rgba8Image{ nullptr }; // src image
    const uint8_t* palette{ nullptr };
    int numColors{ };
    RGB curError; // track error
    RGB lastPixel;   // track last pixel to detect edges

    [[nodiscard]] RGB ClosestPaletteColor(const RGB& color) const
    {
        double best = 1e10;
        RGB ans(0, 0, 0);
        for (int i = 0; i < numColors*3; i += 3)
        {
            RGB t(palette + i);
            const double dist = (t - color).lenSq();
            if (dist < best)
            {
                best = dist;
                ans = t;
            }
        }
        return ans;
    }


    // dither pixel at cur_x , cur_y 
    void DitherPixel()
    {
	    const auto colorPtr = rgba8Image + (xpos + ypos * width) * 4;
	    const RGB pixel(colorPtr);
        constexpr double threshold = 100;
        if ((lastPixel - pixel).lenSq() > threshold)
            curError = RGB();
        lastPixel = pixel;

	    const RGB pixelPlusError = pixel + curError; // the desired current color includes past error
        const RGB closestPaletteColor = ClosestPaletteColor(pixelPlusError); // what we have to actually use
        curError = pixelPlusError - closestPaletteColor; // we need to add error back into picture
        closestPaletteColor.Draw(colorPtr); // draw back into image a valid color
    }

    // direction 0-3 = RDLU
    void DitherAndMove(int direction)
    {
        // dither if in bounds
        if (xpos >= 0 && xpos < width && ypos >= 0 && ypos < height)
            DitherPixel();

        // move to next pixel
        const static int del[] = { 1,0,0,1,-1,0,0,-1 };
        xpos += del[2 * direction];
        ypos += del[2 * direction + 1];
    }

    // direction 0-3 = RDLU
    void HilbertCurve(int level, int direction)
    {
        if (level == 0)
            return;

        const int del = 1 + 2 * (direction & 1);
        const int nxt = (direction + del) % 4;
        const int opp = (direction + 2 * del) % 4;
        const int prv = (direction + 3 * del) % 4;

        HilbertCurve(level - 1, nxt);
        DitherAndMove(opp);
        HilbertCurve(level - 1, direction);
        DitherAndMove(prv);
        HilbertCurve(level - 1, direction);
        DitherAndMove(direction);
        HilbertCurve(level - 1, prv);
    }


public:
    // dither image in place to palette
    void DitherImage(
        int widthIn, int heightIn, // image sizes
        uint8_t * rgba8ImageIn, // src image
        int numColorsIn,
        const uint8_t * paletteIn)
    {

        width = widthIn;
        height = heightIn;
        rgba8Image = rgba8ImageIn;
        palette = paletteIn;
        numColors = numColorsIn;

        /* determine the required order of the Hilbert curve */
        const int size = std::max(width, height);
        const int level = static_cast<int>(std::ceil(std::log2(size)));

        xpos = 0;
        ypos = 0;
        if (level > 0)
            HilbertCurve(level, 3);
        DitherAndMove(0); // dir here doesn't matter
    }
};