#pragma once

#include <cstdlib>

// x (horizontal right), y (vertical down) coordinates, (0, 0) at upper left corner
class BitRaster
{
private:
    typedef unsigned int BYTE_TYPE;
	BYTE_TYPE* bytes;
	int bytesInRow;

	inline int getByteIndex(int x, int y)
	{
		return y * bytesInRow + (x >> BITS_IN_BYTE_POWER_OF_2);
    }

public:
	int w;  // width
	int h;  // height
	static const int BITS_IN_BYTE_POWER_OF_2 = 5;
	static const int BITS_IN_BYTE = 1 << BITS_IN_BYTE_POWER_OF_2;

	BitRaster(int _w, int _h)
	{
		w = _w;
		h = _h;
		bytesInRow = (w - 1) / BITS_IN_BYTE + 1;
		bytes = (BYTE_TYPE*)calloc(bytesInRow * h, sizeof(BYTE_TYPE));
    }
	
	~BitRaster()
	{
		free(bytes);
    }
	
	inline bool getBit(int x, int y)
	{
		return bytes[getByteIndex(x, y)] & (1 << ((BITS_IN_BYTE - 1) & x));
	}

	inline void setBit(int x, int y, bool val)
	{
		if (val) {
			bytes[getByteIndex(x, y)] |= (1 << ((BITS_IN_BYTE - 1) & x));
		}
		else {
			bytes[getByteIndex(x, y)] &= ~(1 << ((BITS_IN_BYTE - 1) & x));
		}
	}

	inline bool byteIsTrue(int x, int y) {
        return bytes[getByteIndex(x, y)] == (1LL << BITS_IN_BYTE) - 1;
	}

	inline bool byteIsFalse(int x, int y) {
		return bytes[getByteIndex(x, y)] == 0;
	}

	inline void setByteTrue(int x, int y)
	{
        bytes[getByteIndex(x, y)] = (1LL << BITS_IN_BYTE) - 1;
    }

	inline void setByteFalse(int x, int y)
	{
		bytes[getByteIndex(x, y)] = 0;
    }
};
