/*
# Copyright(C) 2010-2017 Viacheslav Makoveichuk (email: slv709@gmail.com, skype: vyacheslavm81)
# This file is part of HArray.
#
# HArray is free software : you can redistribute it and / or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HArray is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#ifndef _HARRAY_FIX_BASE		 // Allow use of features specific to Windows XP or later.
#define _HARRAY_FIX_BASE 0x709 // Change this to the appropriate value to target other versions of Windows.
#endif

#define _RELEASE 0x1234567

const uint32_t MAX_CHAR = 256;
const uint32_t MAX_SHORT = 65536;

const uint32_t INIT_MAX_PAGES = 256;
const uint32_t PAGE_SIZE = MAX_SHORT*2;
const uint8_t BLOCK_SIZE = 16;
const uint8_t ROW_LEN = 3;

const uint32_t REPOSITORY_VERSION = 1;

const uint32_t COUNT_TEMPS = 50;

const uint32_t BLOCK_ENGINE_BITS = 4; //bits of block
const uint32_t BLOCK_ENGINE_SIZE = 16; //size of block
const uint32_t BLOCK_ENGINE_SHIFT = 32 - BLOCK_ENGINE_BITS;
const uint32_t BLOCK_ENGINE_STEP = 4;

//const uint32_t BLOCK_ENGINE_BITS = 16; //bits of block
//const uint32_t BLOCK_ENGINE_SIZE = 65536; //size of block
//const uint32_t BLOCK_ENGINE_SHIFT = 32 - BLOCK_ENGINE_BITS;
//const uint32_t BLOCK_ENGINE_STEP = 16;

const uint32_t BRANCH_ENGINE_SIZE = 4; //can be changed

const uint8_t EMPTY_TYPE = 0;
const uint8_t MIN_BRANCH_TYPE1 = 1;
const uint8_t MAX_BRANCH_TYPE1 = BRANCH_ENGINE_SIZE;
const uint8_t MIN_BRANCH_TYPE2 = BRANCH_ENGINE_SIZE + 1;
const uint8_t MAX_BRANCH_TYPE2 = BRANCH_ENGINE_SIZE * 2;
const uint8_t MIN_BLOCK_TYPE = MAX_BRANCH_TYPE2 + 1;
const uint8_t MAX_BLOCK_TYPE = MIN_BLOCK_TYPE + (32 / BLOCK_ENGINE_STEP) - 1;
const uint8_t VAR_TYPE = MAX_BLOCK_TYPE + 1; //for var value
const uint8_t CONTINUE_VAR_TYPE = VAR_TYPE + 1; //for continue var value
const uint8_t CURRENT_VALUE_TYPE = CONTINUE_VAR_TYPE + 1;
const uint8_t VALUE_TYPE = CURRENT_VALUE_TYPE + 1;
const uint8_t ONLY_CONTENT_TYPE = VALUE_TYPE + 1;

const uint8_t MOVES_LEVEL1_STAT = 0;
const uint8_t MOVES_LEVEL2_STAT = 1;
const uint8_t MOVES_LEVEL3_STAT = 2;
const uint8_t MOVES_LEVEL4_STAT = 3;
const uint8_t MOVES_LEVEL5_STAT = 4;
const uint8_t MOVES_LEVEL6_STAT = 5;
const uint8_t MOVES_LEVEL7_STAT = 6;
const uint8_t MOVES_LEVEL8_STAT = 7;
const uint8_t SHORT_WAY_STAT = 8;
const uint8_t LONG_WAY_STAT = 9;
const uint8_t CONTENT_BRANCH_STAT = 10;

const uint8_t CURRENT_VALUE_SEGMENT_TYPE = 1;
const uint8_t BRANCH_SEGMENT_TYPE = CURRENT_VALUE_SEGMENT_TYPE + 1;
const uint8_t BLOCK_VALUE_SEGMENT_TYPE = BRANCH_SEGMENT_TYPE + 1;
const uint8_t BLOCK_BRANCH1_SEGMENT_TYPE = BLOCK_VALUE_SEGMENT_TYPE + 1;
const uint8_t BLOCK_BRANCH2_SEGMENT_TYPE = BLOCK_BRANCH1_SEGMENT_TYPE + 1;
const uint8_t BLOCK_OFFSET_SEGMENT_TYPE = BLOCK_BRANCH2_SEGMENT_TYPE + 1;
const uint8_t VAR_SHUNT_SEGMENT_TYPE = BLOCK_OFFSET_SEGMENT_TYPE + 1;
const uint8_t VAR_VALUE_SEGMENT_TYPE = VAR_SHUNT_SEGMENT_TYPE + 1;

const uint32_t MIN_COUNT_RELEASED_CONTENT_CELLS = MAX_SHORT * 2; //two pages
const uint32_t MIN_COUNT_RELEASED_BRANCH_CELLS = MAX_SHORT;
const uint32_t MIN_COUNT_RELEASED_BLOCK_CELLS = MAX_SHORT;
const uint32_t MIN_COUNT_RELEASED_VAR_CELLS = MAX_SHORT;

const uint8_t MAX_KEY_SEGMENTS = MAX_CHAR - ONLY_CONTENT_TYPE;
const uint8_t MIN_HEADER_BASE_BITS = 14;		  //16384 slots
const uint8_t MAX_HEADER_FILL_FACTOR_BITS = 4; //fill factor 1/16 of header size

typedef bool HARRAY_ITEM_VISIT_FUNC(uint32_t* key, uint32_t keyLen, uint32_t value, void* pData);

struct HArrayFixBaseInfo
{
	uint32_t Version;

	uint32_t KeyLen;
	uint32_t ValueLen;

	uint32_t HeaderBase;

	uint32_t ContentPagesCount;
	uint32_t VarPagesCount;
	uint32_t BranchPagesCount;
	uint32_t BlockPagesCount;

	uint32_t LastContentOffset;
	uint32_t LastVarOffset;
	uint32_t LastBranchOffset;
	uint32_t LastBlockOffset;
};

struct HACursor
{
	uint32_t CountFullContentPage;
	uint32_t SizeLastContentPage;

	uint32_t Page;
	uint32_t Index;

	uint32_t* Value;
};

//struct ContentTypeCell
//{
//	uint8_t Type;
//};

struct HArrayPair
{
public:
	HArrayPair(uint32_t maxKeyLen)
	{
		Key = new uint32_t[maxKeyLen];
		KeyLen = 0;
		Value = 0;
		ValueType = 0;
	}

	uint32_t* Key;
	uint32_t KeyLen;
	uint32_t Value;
	uint8_t ValueType;

	inline void setPair(uint32_t skipKeyLen,
						uint32_t* key,
						uint32_t keyLen,
						uint32_t value,
						uint8_t valueType)
	{
		for (uint32_t i=0; skipKeyLen < keyLen; i++, skipKeyLen++)
		{
			Key[i] = key[skipKeyLen];
		}

		KeyLen = keyLen;
		Value = value;
		ValueType = valueType;
	}

	void print()
	{
		for(uint32_t i = 0; i < KeyLen; i++)
		{
			printf("%u ", Key[i]);
		}

		printf("=> %u\n", Value);
	}

	~HArrayPair()
	{
		delete[] Key;
	}
};

struct BranchCell
{
	uint32_t Values[BRANCH_ENGINE_SIZE];
	uint32_t Offsets[BRANCH_ENGINE_SIZE];
};

struct BlockCell
{
	uint8_t Type;
	uint32_t Offset;
	uint32_t ValueOrOffset;
};

//struct ContentCell
//{
//	uint8_t Type;
//	uint32_t Value;
//};

struct VarCell
{
	uint8_t ValueContCellType;
	uint32_t ValueContCellValue;

	uint8_t ContCellType;
	uint32_t ContCellValue;
};

struct ContentPage
{
	uint8_t pType[MAX_SHORT];
	uint32_t pContent[MAX_SHORT];
};

struct VarPage
{
	VarCell pVar[MAX_SHORT];
};

struct BranchPage
{
	BranchCell pBranch[MAX_SHORT];
};

struct BlockPage
{
	BlockPage()
	{
		memset(this, 0, sizeof(BlockPage));
	}

	BlockCell pBlock[MAX_SHORT];
};

struct CompactPage
{
	CompactPage()
	{
		memset(this, 0, sizeof(CompactPage));
	}

	uint32_t Values[MAX_CHAR];
	uint32_t Offsets[MAX_CHAR];

	uint32_t Count;

	CompactPage* pNextPage;
};

struct SegmentPath
{
	uint8_t Type;

	uint8_t* pContentCellType;
	uint32_t* pContentCellValue;

	BlockCell* pBlockCell;
	uint32_t StartBlockOffset;

	BranchCell* pBranchCell1;
	uint32_t BranchOffset1;

	BranchCell* pBranchCell2;
	uint32_t BranchOffset2;

	VarCell* pVarCell;
	uint32_t VarOffset;

	uint32_t BranchIndex;
	uint32_t ContentOffset;
	uint32_t BlockSubOffset;

	void print()
	{
		/*
		printf("Type: %u, ", Type);

		if(pContentCell)
			printf("ContentCell: Type=%u, Value=%u, ", pContentCell->Type, pContentCell->Value);

		if(pBlockCell)
			printf("BlockCell: Type=%u, Offset=%u, ValueOrOffset=%u, ", pBlockCell->Type, pBlockCell->Offset, pBlockCell->ValueOrOffset);

		printf("StartBlockOffset: %u, ", StartBlockOffset);

		if (pBranchCell1)
			printf("pBranchCell1: Type=%u, Offset=%u, ValueOrOffset=%u, ", pBlockCell->Type, pBlockCell->Offset, pBlockCell->ValueOrOffset);
		*/
	}
};

typedef uint32_t (*NormalizeFunc)(void* key);

typedef int (*CompareFunc)(void* key1, uint32_t keyLen1, void* key2, uint32_t keyLen2);

typedef int(*CompareSegmentFunc)(void* keySeg1, void* keySeg2, uint32_t index);
