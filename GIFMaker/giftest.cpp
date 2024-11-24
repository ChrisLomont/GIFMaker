// Chris Lomont 2022
// simplest GIF writer

#include <vector>
#include <iostream>
#include <filesystem>
#include "GIFMaker.h"

namespace fs = std::filesystem;

using namespace std;

void Fill(vector<uint8_t> & rgba, int w, int h, int pct)
{
	int ww = w * pct / 100, hh = h * pct / 100;
	cout << ww << ' ' << hh << endl;
	for (int j = 0; j < hh; ++j)
		for (int i = 0; i < ww; ++i)
		{
			auto index = (i + j * w) * 4;
#if 0
			rgba[index++] = frame * i + frame * 3 * j; // r
			rgba[index++] = 2*frame * i + frame * 1 * j; // g
			rgba[index++] = frame * i + frame * 2 * j; // b
#else
			rgba[index++] = 255;
			rgba[index++] = 0;
			rgba[index++] = 255*pct/100;
#endif
			rgba[index] = 255; // alpha
		}
}


void DumpB(int bits, vector<uint8_t> dat, int w, int h)
{
	int bitPos = 0, symbs = 0;

	int x = 0, y = 0;


	while (bitPos<dat.size()*8-7)
	{
		uint32_t symb = 0;

		for (int i = 0; i < bits; ++i)
		{
			uint8_t c = dat[bitPos / 8];
			uint8_t mask = 1 << (bitPos & 7);
			uint8_t bit = (c & mask) ? 1 : 0;
			symb |= (bit << i);
			++bitPos;
		}
		cout << symb << " ";
		symbs++;
		if (symb < 255)
		{
			x++; if (x == w) {
				x = 0; cout << endl;
				++y;
				cout << y << "> ";
			}
		}
	}
	cout << endl;
	cout << "Symbols " << symbs << endl;
}

void Decode(const std::string& filename, int w, int h)
{
	// decode LZW chunk

	fs::path inputFilePath{ filename };
	auto length = fs::file_size(inputFilePath);
	if (length == 0) {
		return;
	}
	std::vector<uint8_t> buffer(length);
	std::ifstream inputFile(filename, std::ios_base::binary);
	inputFile.read(reinterpret_cast<char*>(buffer.data()), length);
	inputFile.close();


	// byte 2C starts block
	int count = 0, start = 0, end = 0;
	for (int i = 0; i < buffer.size(); ++i)
	{
		auto b = buffer[i];
		if (b == 0x2C)
		{
			++count;
			start = i;
		}
		if (b == 0x3B && start > 0)
			end = i;
	}
	cout << "2c count " << count << " start:end " << start << ":" << end << endl;

	vector<uint8_t> lzwbits;

	start += 10; // point to LZW part
	int p = start;
	int codeSize = buffer[p++];
	cout << "code size " << codeSize << endl;
	while (true)
	{
		int blockLen = buffer[p++];
		cout << "block length " << blockLen << endl;
		for (int i = 0; i < blockLen; ++i)
			lzwbits.push_back(buffer[p++]);
		if (blockLen == 0) 
			break;
	}

	// now dump 9 bits at a time
	DumpB(9,lzwbits, w, h);

	// 0 terminator, then end?
	cout << "Block term 0 = " << (int)buffer[p++] << endl;
	cout << "pos now " << p << endl;
	cout << "File end 0x3B (59) = " << (int)buffer[p++] << endl;
}

int main()
{
	int s = 100, maxf = 10;
	int w = s, h = s;
	vector<uint8_t> pix(w * h * 4);
	memset(pix.data(),0,pix.size());

	int delay = 50; // 100ths of sec

	GIFWriter gif;
	string filename = "test.gif";
	gif.Start(filename, w, h);
	for (int f = 0; f < maxf; ++f)
	{
		int pct = maxf > 1 ? f * 100 / (maxf - 1) : 50;
		Fill(pix,w,h,pct);
		gif.AddFrame(pix.data(), delay);
	}
	gif.Write();
	cout << "File " << filename << " written\n";
	//Decode(filename, w, h);
	return 0;
}
