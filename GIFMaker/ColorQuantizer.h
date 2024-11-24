#pragma once

// Chris Lomont 2022
// color quantization

#include <algorithm>
#include <numeric>
#include <vector>

class ColorQuantizer {

    struct RGB {
        RGB(double r=0, double g=0, double b=0) :r(r),g(g),b(b){}
        RGB(const uint8_t * rgba, int index) : RGB(rgba[4 * index], rgba[4 * index+1], rgba[4 * index+2]){}
        void Out3(uint8_t* rgba) const { rgba[0] = r; rgba[1] = g; rgba[2] = b; }
    	double r{ 0 }, g{ 0 }, b{ 0 };
    	friend RGB operator*(double s, const RGB& rgb){ return RGB{ s * rgb.r,s * rgb.g, s * rgb.b }; }
        friend RGB operator-(const RGB& lhs, const RGB& rhs) { return RGB{lhs.r - rhs.r,lhs.g - rhs.g,lhs.b - rhs.b};}
        friend RGB operator+(const RGB& lhs, const RGB& rhs) { return RGB{lhs.r + rhs.r,lhs.g + rhs.g,lhs.b + rhs.b};}
        [[nodiscard]] double lengthSq() const { return r * r + g * g + b * b;}
    };

    struct Cluster { RGB center; int size{}; };

	// 2D quasi random Sobel sequence, Numerical Recipes in C
    // pt in [0,1)x[0,1)
    class Sobel2D
    {
        enum { BIT_SIZE = 30 };
        uint64_t ix{}, iy{}, index{0}, * iu[2 * BIT_SIZE + 1]{};
        uint64_t iv[2 * BIT_SIZE + 1] =
        { 0, 1, 1, 1, 1, 1, 1, 3, 1, 3, 3, 1, 1, 5, 7, 7, 3, 3, 5, 15, 11, 5, 15, 13, 9 };
    public:
        Sobel2D()
        {
            constexpr int mdeg[3] = { 0, 1, 2 };
            constexpr int ip[3] = { 0, 0, 1 };
            for (int j = 1, k = 0; j <= BIT_SIZE; j++, k += 2) 
                iu[j] = &iv[k];

            for (int k = 1; k <= 2; k++) {
                for (int j = 1; j <= (int)mdeg[k]; j++) 
                    iu[j][k] <<= (BIT_SIZE - j);

                for (int j = mdeg[k] + 1; j <= BIT_SIZE; j++) {
                    uint64_t ipp = ip[k];
                    uint64_t i = iu[j - mdeg[k]][k];
                    i ^= (i >> mdeg[k]);

                    for (int l = mdeg[k] - 1; l >= 1; l--) {
                        if (ipp & 1) {
                            i ^= iu[j - l][k];
                        }

                        ipp >>= 1;
                    }

                    iu[j][k] = i;
                }
            }
        }
        std::pair<double,double> Next()
        {
            constexpr double scale = 1.0 / (1L << BIT_SIZE);
            uint64_t tmp = index++, i=1;
            while (i <= BIT_SIZE && (tmp&1)) {
                tmp >>= 1;
                ++i;
            }
            const auto k = (i - 1) * 2;
            ix ^= iv[k + 1];
            iy ^= iv[k + 2];
            return { ix * scale,iy * scale };
        }
    };

    int width{ 0 }, height{ 0 }, size{ 0 };
    Sobel2D sobel2D;
    const uint8_t* inPixelsRgba{ nullptr };

    // find index of closest cluster
    static int ClosestIndex(const RGB & rgb, int colorCount, const Cluster* clusters)
    {
        double minDist = std::numeric_limits<double>::max();
        int minIndex = 0; // assume
        for (int j = 0; j < colorCount; j++)
        {
            const double dist = (clusters[j].center - rgb).lengthSq();

            if (dist < minDist)
            {
                minDist = dist;
                minIndex = j;
            }
        }
        return minIndex;
    }

    // "Fast Color Quantization Using MacQueen's K-Means Algorithm"
    // S. Thompson, M. E. Celebi, and K. H. Buck, 2020
    void OnlineKMeans(const int colorCount, Cluster* clusters)
    {
        const auto sampleCount = lround(0.5 * size);

        for (int i = 0; i < sampleCount; i++)
        {
            auto [rx,ry] = sobel2D.Next();

            const int row = static_cast<int>(std::clamp(round(ry * height), 0.0, height - 1.0));
            const int col = static_cast<int>(std::clamp(round(rx * width), 0.0, width - 1.0));

            RGB sampledPixel(inPixelsRgba, row * width + col);

            const int minIndex = ClosestIndex(sampledPixel, colorCount, clusters);

            const int oldSize = clusters[minIndex].size;
            const int newSize = oldSize + 1;

            const double learningRate = pow(newSize, -0.5);

            clusters[minIndex].center = clusters[minIndex].center + learningRate * (sampledPixel - clusters[minIndex].center);
            clusters[minIndex].size = newSize;
        }
    }
    

    // "Incremental Online K-Means Algorithm",
    // A. D. Abernathy and M. E. Celebi, 2022
    void IncrementalOnlineKMeans(std::vector<Cluster>& clusters)
    {
        const int colorCount = clusters.size();
        const int splitCount = static_cast<int>(lround(log2(colorCount)));

        std::vector<Cluster> tempClusters((2 * colorCount - 1));
        
    	RGB sum(0, 0, 0); // centroid
        for (int i = 0; i < size; i++) sum = sum + RGB(inPixelsRgba, i);
        tempClusters[0].center = (1.0 / size) * sum;
        tempClusters[0].size = 0;

        for (int t = 0; t < splitCount; t++)
        {
	        const auto t2 = 1 << t;
            for (int n = t2 - 1; n < 2*t2 - 1; n++)
            {
                const auto pixel = tempClusters[n].center;
                tempClusters[2 * n + 1] = { pixel,0 };
                tempClusters[2 * n + 2] = { pixel,0 };
            }
            OnlineKMeans(2*t2,tempClusters.data() + 2*t2 - 1);
        }
        for (int j = 0; j < colorCount; j++) {
            clusters[j].center = tempClusters[j + colorCount - 1].center;
        }
    }

public:
    /**
     * \brief Color quantizer
     * \tparam T type for output indices, deduced from outPixelIndices 
     * \param numColors number of colors to reduce to, must be power of 2 in 2 to 65536
     * \param widthIn width of input in pixels
     * \param heightIn height of input in pixels
     * \param inPixelsRgbaIn RGBA8 input pixels, row major
     * \param outPixelIndices indices into palette of output
     * \param paletteRgb output palette RGB8 (NOTE: alpha not in here, must be dealth with outside this)
     */
    template<typename T>
    void RemapImage(
        int numColors, // power of 2 in 2 to 65536
        const int widthIn, const int heightIn,
        const uint8_t* inPixelsRgbaIn, // RGBA, row major
        T * outPixelIndices,      // indices into paletteRgb
        uint8_t* paletteRgb
    )
    {
        width = widthIn;
        height = heightIn;
        size = width * height;
        inPixelsRgba = inPixelsRgbaIn;

        sobel2D = Sobel2D(); // make new one

        // generate palette
        std::vector<Cluster> palette(numColors);

        IncrementalOnlineKMeans(palette);

        for (int i = 0; i < size; i++)
            outPixelIndices[i] = T(ClosestIndex(RGB(inPixelsRgba, i), numColors, palette.data()));

        // convert palette
        for (int i = 0; i < numColors; ++i)
            palette[i].center.Out3(paletteRgb+3*i);
    }
};
