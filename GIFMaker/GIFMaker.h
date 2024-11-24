#pragma once
// single file, simple animated GIF file generator
// modern C++
// Chris Lomont 2022

#include <string>
#include <cstdint>
#include <unordered_map>
#include <fstream>
#include <array>
// GIF writer writes GIFs from RGBA byte buffer frames
//
class GIFWriter
{
	// good GIF format description
	// https://giflib.sourceforge.net/whatsinagif/index.html
	// also the gif89a.txt standard, found online

	std::ofstream os;
	int w{ 0 }, h{0};
	int bitDepth = 8; // 8 bits palette - smaller? todo
public:
	int useLocalPalettes = 0;  // 0 or 1 todo
	int useGlobalPalettes = 1; // 0 or 1 todo
	int transparentIndex = 255; // user sets? todo

	// palettes https://fornaxvoid.com/colorpalettes/
	// Deluxe Paint 256 color palette
	// todo - can clamp to mults of 16, pack it more
	std::vector<uint8_t> palette{
0, 0, 0, 48, 138, 69, 227, 227, 227, 113, 113, 113, 243, 0, 0, 146,
0, 0, 243, 211, 211, 243, 162, 81, 243, 243, 211, 227, 211, 0, 211,
243, 81, 211, 243, 211, 0, 243, 0, 0, 146, 0, 211, 243, 243, 0, 227,
227, 170, 170, 170, 69, 223, 69, 223, 223, 223, 113, 113, 113, 227,
0, 0, 130, 0, 0, 243, 178, 178, 243, 146, 65, 243, 243, 178, 195,
195, 0, 195, 243, 65, 178, 243, 178, 0, 227, 0, 0, 130, 0, 178, 243,
243, 0, 195, 195, 101, 101, 101, 69, 223, 207, 223, 223, 223, 101,
101, 101, 227, 0, 0, 113, 0, 0, 243, 146, 146, 243, 130, 32, 243,
243, 146, 178, 162, 0, 178, 243, 32, 146, 243, 146, 0, 227, 0, 0,
113, 0, 146, 243, 243, 0, 178, 178, 223, 223, 223, 48, 138, 207, 195,
195, 195, 81, 81, 81, 211, 0, 0, 113, 0, 0, 243, 113, 113, 243, 113,
0, 243, 243, 113, 146, 146, 0, 162, 243, 0, 130, 243, 113, 0, 211, 0,
0, 113, 0, 113, 243, 243, 0, 146, 146, 207, 48, 69, 138, 138, 223,
178, 178, 178, 65, 65, 65, 195, 0, 0, 97, 0, 0, 243, 81, 81, 227, 97,
0, 243, 243, 81, 130, 130, 0, 146, 227, 0, 97, 243, 81, 0, 195, 0, 0,
97, 0, 81, 243, 243, 0, 130, 130, 223, 138, 69, 69, 48, 207, 170,
170, 170, 48, 48, 48, 178, 0, 0, 81, 0, 0, 243, 65, 65, 195, 97, 0,
243, 243, 65, 113, 97, 0, 130, 195, 0, 65, 243, 65, 0, 178, 0, 0, 81,
0, 65, 243, 243, 0, 113, 113, 207, 223, 69, 207, 48, 207, 146, 146,
146, 32, 32, 32, 178, 0, 0, 65, 0, 0, 243, 32, 32, 178, 81, 0, 243,
243, 32, 81, 81, 0, 113, 178, 0, 32, 243, 32, 0, 178, 0, 0, 65, 0,
32, 243, 243, 0, 81, 81, 138, 138, 48, 223, 138, 207, 130, 130, 130,
32, 32, 32, 162, 0, 0, 65, 0, 0, 243, 0, 0, 146, 65, 0, 243, 243, 0,
65, 65, 0, 97, 146, 0, 0, 243, 0, 0, 162, 0, 0, 65, 0, 0, 243, 243,
0, 65, 65, 81, 178, 243, 211, 211, 243, 0, 0, 243, 0, 0, 146, 243,
211, 243, 146, 0, 227, 243, 211, 243, 227, 0, 227, 243, 227, 211,
195, 146, 130, 130, 65, 48, 81, 16, 16, 243, 81, 81, 195, 32, 32, 32,
195, 48, 81, 32, 195, 65, 178, 243, 178, 178, 243, 0, 0, 227, 0, 0,
130, 227, 178, 243, 130, 0, 195, 243, 178, 243, 195, 0, 195, 243,
211, 211, 178, 130, 113, 130, 48, 48, 65, 16, 0, 243, 178, 130, 195,
65, 32, 32, 195, 81, 130, 32, 195, 32, 162, 243, 146, 146, 243, 0, 0,
227, 0, 0, 113, 211, 146, 243, 113, 0, 178, 243, 146, 243, 178, 0,
178, 243, 211, 195, 178, 113, 97, 113, 48, 32, 65, 0, 0, 243, 243,
130, 195, 113, 32, 32, 195, 130, 178, 32, 195, 0, 146, 243, 113, 130,
243, 0, 0, 211, 0, 0, 113, 211, 113, 243, 97, 0, 146, 243, 113, 243,
146, 0, 146, 227, 195, 178, 162, 113, 97, 113, 48, 32, 48, 0, 0, 130,
243, 130, 195, 146, 32, 32, 195, 178, 195, 32, 162, 0, 130, 227, 81,
97, 243, 0, 0, 195, 0, 0, 97, 195, 81, 243, 81, 0, 130, 243, 81, 243,
130, 0, 130, 227, 178, 162, 162, 97, 81, 113, 32, 32, 48, 0, 0, 130,
243, 243, 195, 195, 32, 32, 162, 195, 195, 32, 130, 0, 113, 195, 65,
65, 243, 0, 0, 178, 0, 0, 81, 178, 65, 243, 65, 0, 113, 243, 65, 243,
97, 0, 113, 211, 178, 146, 162, 97, 81, 97, 32, 16, 32, 0, 0, 130,
130, 243, 146, 195, 32, 32, 113, 195, 195, 32, 81, 0, 97, 178, 32,
32, 243, 0, 0, 178, 0, 0, 65, 178, 32, 243, 48, 0, 81, 243, 32, 243,
81, 0, 81, 211, 162, 146, 146, 81, 65, 97, 16, 16, 32, 0, 0, 178,
130, 243, 113, 195, 32, 32, 81, 195, 195, 32, 32, 0, 81, 146, 0, 0,
243, 0, 0, 162, 0, 0, 65, 162, 0, 243, 32, 0, 65, 243, 0, 243, 65, 0,
65, 195, 146, 130, 146, 81, 65, 81, 16, 16, 32, 0, 0, 243, 130, 243,
65, 195, 32, 32, 32, 195, 243, 243, 243 };


	bool Start(const std::string& filename, int w, int h)
	{
		os.open(filename, std::ios::out | std::ios::binary);
		if (!os.is_open())
			return false;
		this->w = w;
		this->h = h;

		// Header Block
		Write("GIF89a");

		// Logical Screen Descriptor
		Write16(w); // canvas width
		Write16(h); // canvas height
		Write8(
			(useGlobalPalettes << 7) | // 1: has global color tbl 
			(((bitDepth-1)&7) << 4)  | // 3: color res this + 1 bits / channel
			(0 << 3)                 | // 1: palette not sorted
			((bitDepth-1)&7)           // 3: table size 2^(this+1)
		);
		// background color index, pixel aspect ratio not mentioned (assume 1.0)
		Write8({ 0,0 });  

		// Global Color Table
		if (useGlobalPalettes) // todo - let user set?
		{
			os.write((const char*)palette.data(),palette.size());
		}

		// Application Extension: animation header (todo - only if asked?)
		// GIF Extension code 0x21,  App specific 0xFF, length 0x0B
		Write8({0x21,0xFF,0x0B});
		Write("NETSCAPE2.0"); // legacy :)
		// 3 bytes data, must be 1, times to loop Lo/Hi, block terminator
		Write8({3,1,0,0,0});

		// Comment block to tag this header
		std::string comment = "LomontGIF1.0";

		// GIF Extension code 0x21, Comment 0xFE, txt length
		Write8({ 0x21,0xFE,(int)(comment.size())});
		Write(comment.c_str()); // comment text
		Write8(0);    // block terminator
	}


	// delay in 100ths of sec
	// frame is taken so first 255 unique colors are palette, rest mapped to those
	void AddFrame(const uint8_t * rgbaData, uint16_t delayIn100thsOfSecond)
	{
		// todo - track alpha, local color table
		std::vector<uint8_t> pixels;
		ColorMapper(rgbaData, pixels);


		int pack =
			(0 << 5) | // 3: reserved
			(1 << 2) | // 3: disposal method leave in place
			(0 << 1) | // 1: no user input
			(0 << 0);   // 1: no transparency
		int delayLo = delayIn100thsOfSecond & 255;
		int delayHi = delayIn100thsOfSecond/256;
		// Extension 0x21, Graphic Control Extension 0xF9, byte size 4,
		// delay lo byte, delay high byte, transparent color index, block end
		Write8({ 0x21,0xF9,0x04, pack, delayLo, delayHi, transparentIndex, 0});


		// Image Descriptor 0x2C, left lo, left hi, top lo, top hi,
		// width low, width hi, height lo, height hi
		Write8({ 0x2C,0,0,0,0,w & 255,w / 256,h & 255,h / 256 });
		Write8(
			((useLocalPalettes ?1:0) << 7) | // 1: has local color table
			(0 << 6) |                       // 1: interlaced
			(0 << 5) |                       // 1: sorted palette
			(0 << 3) |                       // 2: reserved
			((bitDepth-1)&7)                 // 3: table size
		);
		if (useLocalPalettes)
		{
#if 0
			// Local color table
			index = 0;
			while (index < 256)
			{
				uint32_t c = index < colors.size() ? colors[index] : 0;
				Write8(c >> 24); // r
				Write8(c >> 16); // g
				Write8(c >> 8); // b
				++index;
			}
#endif
		}


		// LZW table data
		LZW lzw(os,8);
		for (auto p : pixels)
		{
			lzw.AddPixelIndex(p);
		}
		lzw.End();

		Write8(0);   // block end

	}
	bool Write()
	{
		// trailer marks end of file
		Write8(0x3B);
		os.close();
		return true;
	}

private:
	// color quantization survey
	// "Forty years of color quantization: a modern, algorithmic survey," 2023, M. Emre Celebi
	void WritePalette()
	{
//		todo
	}

	// map data into pixel indices
	void ColorMapper(const uint8_t* rgbaData, std::vector<uint8_t> & pixels)
	{
		pixels.clear();
		int index = 0;
		while (index < w * h * 4)
		{
			uint32_t r = rgbaData[index++];
			uint32_t g = rgbaData[index++];
			uint32_t b = rgbaData[index++];
			uint32_t a = rgbaData[index++];

			// closest color
			pixels.push_back(GetIndex(r, g, b, a));
		}		
	}

	void Write(const char* str)
	{
		while (*str != 0)
			Write8(*str++);
	}
	void Write8(int v)
	{
		os.put(v);
	}
	void Write8(std::initializer_list<int> lst)
	{
		for (auto v : lst)
			os.put(v);
	}

	void Write16(int v)
	{
		os.put(v & 255);
		os.put((v >> 8));
	}
	void Write(int val, int count)
	{
		for (int i = 0; i < count; ++i)
			Write8(val);
	}

	// can fake out table so it's not compressed,
	// https://github.com/Distrotech/libungif/blob/master/UNCOMPRESSED_GIF
	// but decoders still work
	// idea is to write each byte index as a 9 bit value, then every so often, emit clear code
	struct LZW
	{
		std::ostream& os;
		int codeSize=8;
		LZW(std::ostream& os, int codeSize) : os(os), codeSize(codeSize)
		{
			os.put(codeSize); // min code size
			ClearCode();
		}

		void ClearCode()
		{
			WriteBits(1 << (codeSize)); // clear code 0b1_0000_0000
		}
		void EndCode()
		{
			WriteBits((1 << (codeSize))+1); // end code 0b1_0000_0001
		}


		// data here, each block at most 255 bytes in length
		std::array<uint8_t,255> dat{0};

		uint32_t bitPos = 0; // next open bit position (skip initial byte)

		// bit writing:
		void WriteBit(int val)
		{ // fill in bits low to high, flush, continue
			uint32_t shift = bitPos & 7;
			uint32_t bit = (val & 1) << shift; 
			uint32_t byteIndex = bitPos >> 3;
			dat.at(byteIndex) |= bit;
			bitPos++;
			if (bitPos == 8*255) {
				FlushBlock(); }
		}

		void WriteBits(int val)
		{
			for (int i = 0; i < codeSize+1; ++i)			
				WriteBit((val>>i) & 1);			
		}

		void FlushBlock()
		{
			const int byteCount = bitPos / 8;
			os.put(byteCount); // length
			os.write((const char*)dat.data(), byteCount); // block
			bitPos = 0; // reset
			std::memset(dat.data(), 0, dat.size());
		}

		int entries = 0; // tbl entries since last clear
		// add the next pixel index to the LZW
		void AddPixelIndex(int pixelIndex)
		{
			WriteBits(pixelIndex);
			entries++;
			// todo - this too often, must deduce the rule?
			// 256 initial symbols 0-255, plus end,clear = 258
			// when hits 512, gets bit, so cannot hit that
			if (entries > 200)
			{
				entries = 0;
				ClearCode(); 
			}
		}

		// end the LZW image
		void End()
		{
			ClearCode(); // todo - is needed?			
			EndCode(); // end of frame

			// round up to pad byte, then flush
			bitPos = ((bitPos + 7) / 8) * 8;			
			FlushBlock();
			os.put(0); // 0 length block
		}		
	};


	// todo - make all this better, fix LRU code
	std::vector<uint32_t> lruColors;
	// slow :)
	int GetIndex(int r, int g, int b, int a)
	{
		// todo - add alpha stuff
		uint32_t packed = (r << 24) | (g << 16) | (b << 8); // ignore alpha
		// check lru:
		for (int i = 0; i < lruColors.size(); ++i)
		for (uint32_t v : lruColors)
		{
			uint32_t c = v & 0xFFFFFF00U;
			if (c == packed)
			{
				// matches
				return v & 255; 
			}
		}

		// walk all
		int bestIndex = -1;
		int bestScore = 256 * 256 * 4; // way bigger
		for (int p = 0; p < 256; ++p)
		{
			auto dr = palette[p * 3] - r;
			auto dg = palette[p * 3 + 1] - g;
			auto db = palette[p * 3 + 2] - b;
			auto d = dr * dr + dg * dg + db * db;
			if (d < bestScore)
			{
				bestScore = d;
				bestIndex = p;
			}
		}
		if (lruColors.size() < 10)
		{
			lruColors.push_back(packed | bestIndex);			
		}
		return bestIndex;
	}


#if 0
	/* Color quantization using Macqueen's k-means algorithm */
/*
  For detailed information, see
  S. Thompson, M. E. Celebi, and K. H. Buck,
  Fast Color Quantization Using Macqueens K-Means Algorithm,
  Journal of Real-Time Image Processing, to appear
  (https://doi.org/10.1007/s11554-019-00914-6), 2020.
 */
	typedef struct
	{
		double red, green, blue;
	} RGB_Pixel;

	typedef struct
	{
		int size;
		RGB_Pixel center;
	} RGB_Cluster;
	typedef struct
	{
		int width, height;
		int size;
		RGB_Pixel* data;
	} RGB_Image;
	typedef unsigned char uchar;
	typedef unsigned long ulong;


	ulong
		genrand_int32(void) { todo; }
	

	/* Function for generating a bounded random integer between 0 and RANGE */
	/* Source: http://www.pcg-random.org/posts/bounded-rands.html */

	uint32_t
		bounded_rand(const uint32_t range)
	{
		uint32_t x = genrand_int32();
		uint64_t m = ((uint64_t)x) * ((uint64_t)range);
		uint32_t l = (uint32_t)m;

		if (l < range)
		{
			uint32_t t = -range;

			if (t >= range)
			{
				t -= range;
				if (t >= range)
				{
					t %= range;
				}
			}

			while (l < t)
			{
				x = genrand_int32();
				m = ((uint64_t)x) * ((uint64_t)range);
				l = (uint32_t)m;
			}
		}

		return m >> 32;
	}
#define MAXBIT 30
	/*
	  Returns two quasirandom numbers from a 2D Sobol
	  sequence. Adapted from Numerical Recipies in C.
	 */

	void
		sob_seq(double* x, double* y)
	{
		int j, k, l;
		ulong i, im, ipp;

		/*
		  The following variables are static since we want their
		  values to remain stored after the function returns. These
		  values represent the state of the quasirandom number generator.
		 */

		static double fac;
		static int init = 0;
		static ulong ix1, ix2;
		static ulong in, * iu[2 * MAXBIT + 1];
		static ulong mdeg[3] = { 0, 1, 2 };
		static ulong ip[3] = { 0, 0, 1 };
		static ulong iv[2 * MAXBIT + 1] =
		{ 0, 1, 1, 1, 1, 1, 1, 3, 1, 3, 3, 1, 1, 5, 7, 7, 3, 3, 5, 15, 11, 5, 15, 13, 9 };

		/* Initialize the generator the first time the function is called */
		if (!init)
		{
			init = 1;
			for (j = 1, k = 0; j <= MAXBIT; j++, k += 2)
			{
				iu[j] = &iv[k];
			}

			for (k = 1; k <= 2; k++)
			{
				for (j = 1; j <= mdeg[k]; j++)
				{
					iu[j][k] <<= (MAXBIT - j);
				}

				for (j = mdeg[k] + 1; j <= MAXBIT; j++)
				{
					ipp = ip[k];
					i = iu[j - mdeg[k]][k];
					i ^= (i >> mdeg[k]);

					for (l = mdeg[k] - 1; l >= 1; l--)
					{
						if (ipp & 1)
						{
							i ^= iu[j - l][k];
						}

						ipp >>= 1;
					}

					iu[j][k] = i;
				}
			}

			fac = 1.0 / (1L << MAXBIT);
			in = 0;
		}

		/* Now calculate the next pair of numbers in the 2-D Sobol sequence */

		im = in;
		for (j = 1; j <= MAXBIT; j++)
		{
			if (!(im & 1))
			{
				break;
			}

			im >>= 1;
		}

		im = (j - 1) * 2;
		*x = (ix1 ^= iv[im + 1]) * fac;
		*y = (ix2 ^= iv[im + 2]) * fac;

		in++;

		/* X and Y will fall in [0,1] */
	}

	/* Maximin initialization method */
/*
   For a comprehensive survey of k-means initialization methods, see
   M. E. Celebi, H. Kingravi, and P. A. Vela,
   A Comparative Study of Efficient Initialization Methods for the K-Means Clustering Algorithm,
   Expert Systems with Applications, vol. 40, no. 1, pp. 200210, 2013.
 */
 /* Maximum possible RGB distance = 3 * 255 * 255 */
#define MAX_RGB_DIST 195075 
	void
		maximin(const RGB_Image* img, RGB_Cluster* clusters, const int num_colors, const RGB_Pixel* mean)
	{
		int i, j, max_dist_index = 0;
		double delta_red, delta_green, delta_blue;
		double dist, max_dist;
		double* nc_dist;
		RGB_Pixel pixel;
		RGB_Cluster* cluster;

		nc_dist = (double*)malloc(img->size * sizeof(double));

		/* Initialize first center by the mean R, G, B color */
		cluster = &clusters[0];
		cluster->center.red = mean->red;
		cluster->center.green = mean->green;
		cluster->center.blue = mean->blue;
		cluster->size = 1;

		/* Initialize the nearest-center-distance for each pixel */
		for (j = 0; j < img->size; j++)
		{
			nc_dist[j] = MAX_RGB_DIST;
		}

		/* Choose the remaining centers using maximin */
		for (i = 0 + 1; i < num_colors; i++)
		{
			max_dist = -MAX_RGB_DIST;

			for (j = 0; j < img->size; j++)
			{
				/* Cache the pixel */
				pixel = img->data[j];

				/* Compute the pixel's distance to the previously chosen center */
				cluster = &clusters[i - 1];
				delta_red = pixel.red - cluster->center.red;
				delta_green = pixel.green - cluster->center.green;
				delta_blue = pixel.blue - cluster->center.blue;
				dist = delta_red * delta_red + delta_green * delta_green + delta_blue * delta_blue;

				if (dist < nc_dist[j])
				{
					/* Update the nearest-center-distance for this pixel */
					nc_dist[j] = dist;
				}

				if (max_dist < nc_dist[j])
				{
					/* Update the maximum nearest-center-distance so far */
					max_dist = nc_dist[j];
					max_dist_index = j;
				}
			}

			/* Pixel with maximum distance to its nearest center is chosen as a center */
			cluster = &clusters[i];
			pixel = img->data[max_dist_index];
			cluster->center.red = pixel.red;
			cluster->center.green = pixel.green;
			cluster->center.blue = pixel.blue;
			cluster->size = 1;
		}

		free(nc_dist);
	}

	RGB_Image*
		macqueen_cluster(const RGB_Image* in_img, const int num_colors, const int pres_order,
			const double lr_exp, const double sample_rate, RGB_Pixel* mean)
	{
		int i, j;
		int max_pres, min_dist_index;
		int row_index, col_index, rand_index;
		int old_size, new_size;
		double min_dist, dist;
		double delta_red, delta_green, delta_blue;
		double sob_x, sob_y;
		double rate;
		RGB_Cluster* clusters, * cluster;
		RGB_Pixel in_pix, * out_pix;
		RGB_Image* out_img;

		clusters = (RGB_Cluster*)malloc(num_colors * sizeof(RGB_Cluster));

		auto start = std::chrono::high_resolution_clock::now();

		/* Initialize cluster centers */
		maximin(in_img, clusters, num_colors, mean);

		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop - start);

#ifdef PRINT_TIME_INIT
		printf("Initialization time = %g\n", duration.count() / 1e3);
#endif

		start = std::chrono::high_resolution_clock::now();

		/* Clustering pixel data using Macqueen's k-means algorithm */
		max_pres = in_img->size * sample_rate;
		for (i = 0; i < max_pres; i++)
		{
			/* Choose a pixel quasi- or pseudo-randomly */
			if (pres_order == 0)
			{
				/* Quasirandom */
				sob_seq(&sob_x, &sob_y);

				row_index = (int)(sob_y * in_img->height + 0.5); /* round */
				if (row_index == in_img->height)
				{
					row_index--;
				}

				col_index = (int)(sob_x * in_img->width + 0.5); /* round */
				if (col_index == in_img->width)
				{
					col_index--;
				}

				rand_index = row_index * in_img->width + col_index;
			}
			else
			{
				/* Pseudorandom */
				/* rand_index = ( int ) ( genrand_real2 ( ) * in_img->size ); */
				rand_index = bounded_rand(in_img->size);
			}

			/* Cache the chosen pixel */
			in_pix = in_img->data[rand_index];

			/* Find the nearest center */
			min_dist = MAX_RGB_DIST;
			min_dist_index = -INT_MAX;
			for (j = 0; j < num_colors; j++)
			{
				cluster = &clusters[j];
				delta_red = in_pix.red - cluster->center.red;
				delta_green = in_pix.green - cluster->center.green;
				delta_blue = in_pix.blue - cluster->center.blue;
				dist = delta_red * delta_red + delta_green * delta_green + delta_blue * delta_blue;

				if (dist < min_dist)
				{
					min_dist = dist;
					min_dist_index = j;
				}
			}

			/* Update the nearest center */
			cluster = &clusters[min_dist_index];
			old_size = cluster->size;
			new_size = old_size + 1;
			rate = pow(new_size, -lr_exp);
			cluster->center.red += rate * (in_pix.red - cluster->center.red);
			cluster->center.green += rate * (in_pix.green - cluster->center.green);
			cluster->center.blue += rate * (in_pix.blue - cluster->center.blue);
			cluster->size = new_size;
		}

		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::microseconds> (stop - start);
#ifdef PRINT_TIME_CLUST
		printf("Clustering time = %g\n", duration.count() / 1e3);
#endif

		out_img = (RGB_Image*)malloc(sizeof(RGB_Image));
		out_img->data = (RGB_Pixel*)malloc(in_img->size * sizeof(RGB_Pixel));
		out_img->width = in_img->width;
		out_img->height = in_img->height;
		out_img->size = in_img->size;

		start = std::chrono::high_resolution_clock::now();

		/* Now quantize the image */
		for (i = 0; i < in_img->size; i++)
		{
			/* Cache the pixel */
			in_pix = in_img->data[i];

			/* Find the nearest center */
			min_dist = MAX_RGB_DIST;
			min_dist_index = -INT_MAX;
			for (j = 0; j < num_colors; j++)
			{
				cluster = &clusters[j];

				delta_red = in_pix.red - cluster->center.red;
				delta_green = in_pix.green - cluster->center.green;
				delta_blue = in_pix.blue - cluster->center.blue;
				dist = delta_red * delta_red + delta_green * delta_green + delta_blue * delta_blue;

				if (dist < min_dist)
				{
					min_dist = dist;
					min_dist_index = j;
				}
			}

			/* Replace the input color with the nearest color in the palette */
			in_pix = clusters[min_dist_index].center;
			out_pix = &out_img->data[i];
			out_pix->red = in_pix.red;
			out_pix->green = in_pix.green;
			out_pix->blue = in_pix.blue;
		}

		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::microseconds> (stop - start);
#ifdef PRINT_TIME_MAP
		printf("Mapping time = %g\n", duration.count() / 1e3);
#endif

		free(clusters);

		return out_img;
	}
#endif

};
