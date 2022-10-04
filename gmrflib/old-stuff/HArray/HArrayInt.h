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

#include "HArrayBase.h"
struct HeaderCellInt
{
	uint32_t Offset;
	uint16_t Code;
	uint8_t Type;
};

struct DoubleValueCellInt
{
	uint32_t Value1;
	uint32_t Value2;
	uint16_t Code2;
};

struct MultiplyValueCellInt
{
	uint32_t Value;
	uint16_t Code;
};

class DoublePageInt
{
public:
	DoublePageInt()
	{
		memset(this, 0, sizeof(DoublePageInt));
	}

	DoubleValueCellInt pValues[PAGE_SIZE];

	~DoublePageInt()
	{
	}
};

class MultiplyPageInt
{
public:
	MultiplyPageInt()
	{
		memset(this, 0, sizeof(MultiplyPageInt));
	}

	MultiplyValueCellInt pValues[PAGE_SIZE];

	~MultiplyPageInt()
	{
	}
};

class HArrayInt
{
public:
	HArrayInt(uint32_t headerBase = 20)
	{
		memset(this, 0, sizeof(HArrayInt));

		init(headerBase);
	}

	uint32_t HEADER_BASE;
	uint32_t CONTENT_BASE;

	uint32_t HEADER_SIZE;
	uint32_t BLOCK_SIZE;

	uint32_t BLOCK_BITS;

	HeaderCellInt* pHeader;

	DoublePageInt** pDoublePages;
	MultiplyPageInt** pMultiplyPages;

	uint32_t CountDoublePage;
	uint32_t CountMultiplyPage;

	uint32_t SizeDoublePage;
	uint32_t SizeMultiplyPage;

	uint32_t maxDoubleOffset;
	uint32_t maxMultiplyOffset;

	uint16_t lastHeaderCode;

	uint32_t* pIndexes;
	uint32_t* pValues;

	uint32_t* pReleaseCells;
	uint32_t currReleaseCell;

	void init(uint32_t headerBase)
	{
		destroy();

		//clear pointers
		pHeader = 0;
		pDoublePages = 0;
		pMultiplyPages = 0;

		try
		{
			HEADER_BASE = headerBase;
			CONTENT_BASE = 32 - headerBase;

			uint64_t maxKey = 1;
			maxKey <<= ((uint64_t)headerBase + CONTENT_BASE);

			BLOCK_BITS = CONTENT_BASE;
			HEADER_SIZE = (uint32_t)(maxKey >> BLOCK_BITS);
			BLOCK_SIZE = (uint32_t)((maxKey >> headerBase) - 1);

			pHeader = new HeaderCellInt[(uint64_t)HEADER_SIZE + 1];

			SizeDoublePage = INIT_MAX_PAGES;
			SizeMultiplyPage = INIT_MAX_PAGES;

			pDoublePages = new DoublePageInt * [SizeDoublePage];
			pMultiplyPages = new MultiplyPageInt * [SizeMultiplyPage];

			for (uint32_t i = 0; i <= HEADER_SIZE; i++)
			{
				pHeader[i].Type = 0;
			}

			for (uint32_t i = 0; i < INIT_MAX_PAGES; i++)
			{
				pDoublePages[i] = 0;
				pMultiplyPages[i] = 0;
			}

			pDoublePages[0] = new DoublePageInt();
			pMultiplyPages[0] = new MultiplyPageInt();

			pIndexes = new uint32_t[MAX_SHORT];
			pValues = new uint32_t[MAX_SHORT];
			pReleaseCells = new uint32_t[MAX_SHORT];

			CountDoublePage = 1;
			CountMultiplyPage = 1;

			maxDoubleOffset = 0;
			maxMultiplyOffset = 0;

			lastHeaderCode = 1;

			currReleaseCell = 0;

		}
		catch (...)
		{
			destroy();

			throw;
		}
	}

	void resizeDoublePages()
	{
		DoublePageInt** pNewDoublePages = new DoublePageInt* [(uint64_t)SizeDoublePage * 2];

		uint32_t i = 0;

		for (; i < SizeDoublePage; i++)
		{
			pNewDoublePages[i] = pDoublePages[i];
		}

		SizeDoublePage *= 2;

		for (; i < SizeDoublePage; i++)
		{
			pNewDoublePages[i] = 0;
		}

		delete[] pDoublePages;

		pDoublePages = pNewDoublePages;
	}

	void resizeMultiplyPages()
	{
		MultiplyPageInt** pNewMultiplyPages = new MultiplyPageInt* [(uint64_t)SizeMultiplyPage * 2];

		uint32_t i = 0;

		for (; i < SizeMultiplyPage; i++)
		{
			pNewMultiplyPages[i] = pMultiplyPages[i];
		}

		SizeMultiplyPage *= 2;

		for (; i < SizeMultiplyPage; i++)
		{
			pNewMultiplyPages[i] = 0;
		}

		delete[] pMultiplyPages;

		pMultiplyPages = pNewMultiplyPages;
	}

	bool insert(uint32_t key, uint32_t value)
	{
		try
		{
			uint32_t leftPart = (key >> BLOCK_BITS);
			uint32_t rightPart = (key & BLOCK_SIZE);

			HeaderCellInt& headerCell = pHeader[leftPart];

			if (headerCell.Type)
			{
				//Existing block
				uint32_t countValues = 0;

				if (headerCell.Type == 1) //with one element in block
				{
					if (headerCell.Code == rightPart) //update
					{
						headerCell.Offset = value;
						return false;
					}

					//Allocate block
					if (currReleaseCell > 0)
					{
						uint32_t blockOffset = pReleaseCells[--currReleaseCell];
						DoublePageInt* pPage = pDoublePages[blockOffset >> 16];
						DoubleValueCellInt& valueCell = pPage->pValues[blockOffset & 0xFFFF];

						valueCell.Value1 = headerCell.Offset;

						valueCell.Code2 = rightPart;
						valueCell.Value2 = value;

						headerCell.Offset = blockOffset;
						headerCell.Type = 2;

						return true;
					}

					uint32_t currPage = (maxDoubleOffset >> 16);
					DoublePageInt* pPage = pDoublePages[currPage];

					if (!pPage)
					{
						pPage = new DoublePageInt();
						pDoublePages[currPage] = pPage;
						CountDoublePage++;

						if (CountDoublePage == SizeDoublePage)
						{
							resizeDoublePages();
						}
					}

					DoubleValueCellInt& valueCell = pPage->pValues[maxDoubleOffset & 0xFFFF];

					valueCell.Value1 = headerCell.Offset;

					valueCell.Code2 = rightPart;
					valueCell.Value2 = value;

					headerCell.Offset = maxDoubleOffset++;
					headerCell.Type = 2;

					return true;
				}
				else if (headerCell.Type == 2) //with two elements in block
				{
					uint32_t currPage = (headerCell.Offset >> 16);
					DoublePageInt* pPage = pDoublePages[currPage];

					DoubleValueCellInt& valueCell = pPage->pValues[headerCell.Offset & 0xFFFF];

					if (headerCell.Code == rightPart)
					{
						valueCell.Value1 = value;
						return false;
					}

					if (valueCell.Code2 == rightPart)
					{
						valueCell.Value2 = value;
						return false;
					}

					pIndexes[0] = headerCell.Code;
					pValues[0] = valueCell.Value1;

					pIndexes[1] = valueCell.Code2;
					pValues[1] = valueCell.Value2;

					pIndexes[2] = rightPart;
					pValues[2] = value;

					countValues = 3;

					if (currReleaseCell < MAX_SHORT)
					{
						pReleaseCells[currReleaseCell++] = headerCell.Offset;
					}
				}
				else if (headerCell.Type == 3) //with multiply elements in block
				{
					uint32_t offset = headerCell.Offset + rightPart;

					uint32_t currPage = (offset >> 16);
					MultiplyPageInt* pPage = pMultiplyPages[currPage];

					if (!pPage)
					{
						pPage = new MultiplyPageInt();
						pMultiplyPages[currPage] = pPage;
						CountMultiplyPage++;

						if (CountMultiplyPage == SizeMultiplyPage)
						{
							resizeMultiplyPages();
						}

						MultiplyValueCellInt& valueCell = pPage->pValues[offset & 0xFFFF];
						valueCell.Value = value; //insert
						valueCell.Code = headerCell.Code;

						return true;
					}

					MultiplyValueCellInt& valueCell = pPage->pValues[offset & 0xFFFF];

					if (valueCell.Code == 0)
					{
						valueCell.Value = value; //insert
						valueCell.Code = headerCell.Code;
						return true;
					}
					else if (valueCell.Code == headerCell.Code)
					{
						valueCell.Value = value; //update
						return false;
					}

					//Extract block
					pIndexes[0] = rightPart;
					pValues[0] = value;

					countValues = 1;

					//with calculating MultiplyPageInt
					for (uint32_t i = 0; i <= BLOCK_SIZE; i++)
					{
						uint32_t offset = headerCell.Offset + i;

						uint32_t currPage = (offset >> 16);
						MultiplyPageInt* pPage = pMultiplyPages[currPage];

						if (pPage)
						{
							MultiplyValueCellInt& valueCell = pPage->pValues[offset & 0xFFFF];

							if (valueCell.Code == headerCell.Code)
							{
								pIndexes[countValues] = i;
								pValues[countValues] = valueCell.Value;

								countValues++;

								valueCell.Value = 0;
								valueCell.Code = 0;
							}
						}
					}
				}

				lastHeaderCode++;

				if (lastHeaderCode == 0)
					lastHeaderCode = 1;

				//Find offset for block
				uint32_t blockOffset = maxMultiplyOffset;
				/*
				if(blockOffset > BLOCK_SIZE)
				{
					blockOffset -= BLOCK_SIZE;
				}
				*/

				while (true)
				{
					bool isFound = true;

					for (uint32_t i = 0; i < countValues; i++)
					{
						uint32_t offset = blockOffset + pIndexes[i];

						uint32_t currPage = (offset >> 16);
						MultiplyPageInt* pPage = pMultiplyPages[currPage];

						if (!pPage)
						{
							pPage = new MultiplyPageInt();
							pMultiplyPages[currPage] = pPage;
							CountMultiplyPage++;

							if (CountMultiplyPage == SizeMultiplyPage)
							{
								resizeMultiplyPages();
							}
						}

						MultiplyValueCellInt& valueCell = pPage->pValues[offset & 0xFFFF];

						if (valueCell.Code > 0)
						{
							isFound = false;
							break;
						}
					}

					if (isFound)
					{
						//Allocate block
						for (uint32_t i = 0; i < countValues; i++)
						{
							uint32_t offset = blockOffset + pIndexes[i];

							uint32_t currPage = (offset >> 16);
							MultiplyPageInt* pPage = pMultiplyPages[currPage];
							MultiplyValueCellInt& valueCell = pPage->pValues[offset & 0xFFFF];

							valueCell.Value = pValues[i]; //insert
							valueCell.Code = lastHeaderCode;
						}

						headerCell.Offset = blockOffset;
						headerCell.Code = lastHeaderCode;
						headerCell.Type = 3;

						if (blockOffset > maxMultiplyOffset)
						{
							maxMultiplyOffset = blockOffset;
						}

						break;
					}
					else
					{
						blockOffset++;
					}
				}
			}
			else
			{
				headerCell.Offset = value;
				headerCell.Code = rightPart;
				headerCell.Type = 1;
			}

			return true;
		}
		catch (...)
		{
			destroy();

			throw;
		}
	}

	uint32_t getValueByKey(uint32_t key)
	{
		HeaderCellInt& headerCell = pHeader[key >> BLOCK_BITS];
		uint32_t rightPart = (key & BLOCK_SIZE);

		switch (headerCell.Type)
		{
		case 1:
		{
			if (headerCell.Code == rightPart)
			{
				return headerCell.Offset;
			}

			break;
		}
		case 2:
		{
			DoublePageInt* pPage = pDoublePages[headerCell.Offset >> 16];
			if (pPage)
			{
				DoubleValueCellInt& valueCell = pPage->pValues[headerCell.Offset & 0xFFFF];

				if (headerCell.Code == rightPart)
				{
					return valueCell.Value1;
				}
				else if (valueCell.Code2 == rightPart)
				{
					return valueCell.Value2;
				}
			}

			break;
		}
		case 3:
		{
			uint32_t offset = headerCell.Offset + rightPart;

			MultiplyPageInt* pPage = pMultiplyPages[offset >> 16];
			if (pPage)
			{
				MultiplyValueCellInt& valueCell = pPage->pValues[offset & 0xFFFF];

				if (valueCell.Code == headerCell.Code)
				{
					return valueCell.Value;
				}
			}

			break;
		}
		};

		return 0;
	}


	uint32_t getValuesByRange(uint32_t* buffer,
		uint32_t size,
		uint32_t startKey,
		uint32_t endKey)
	{
		uint32_t count = 0;

		uint32_t startLeftPart = (startKey >> BLOCK_BITS);
		uint32_t endLeftPart = (endKey >> BLOCK_BITS);

		//start
		for (uint32_t currHeader = startLeftPart; currHeader <= endLeftPart; currHeader++)
		{
			HeaderCellInt& headerCell = pHeader[currHeader];
			if (headerCell.Type == 1)
			{
				buffer[count++] = headerCell.Offset;
			}
			else if (headerCell.Type == 2)
			{
				DoublePageInt* pPage = pDoublePages[headerCell.Offset >> 16];

				if (pPage)
				{
					DoubleValueCellInt& valueCell = pPage->pValues[headerCell.Offset & 0xFFFF];
					if (headerCell.Code < valueCell.Code2)
					{
						buffer[count++] = valueCell.Value1;
						buffer[count++] = valueCell.Value2;
					}
					else
					{
						buffer[count++] = valueCell.Value2;
						buffer[count++] = valueCell.Value1;
					}
				}
			}
			else if (headerCell.Type == 3)
			{
				uint32_t startBlockOffset = headerCell.Offset;

				if (currHeader == startLeftPart)
				{
					startBlockOffset += (startKey & BLOCK_SIZE);
				}

				uint32_t endBlockOffset = headerCell.Offset;

				if (currHeader == endLeftPart)
				{
					endBlockOffset += (endKey & BLOCK_SIZE);
				}
				else
				{
					endBlockOffset += BLOCK_SIZE;
				}

				for (uint32_t offset = startBlockOffset; offset <= endBlockOffset; offset++)
				{
					uint32_t currPage = (offset >> 16);
					MultiplyPageInt* pPage = pMultiplyPages[currPage];

					if (pPage)
					{
						MultiplyValueCellInt& valueCell = pPage->pValues[offset & 0xFFFF];

						if (valueCell.Code == headerCell.Code)
						{
							buffer[count++] = valueCell.Value;
						}

						if (count >= size)
						{
							break;
						}
					}
				}
			}
		}

		return count;
	}

	void rebuild()
	{
		//create new
		HArrayInt* pNewArray = new HArrayInt();
		pNewArray->init(HEADER_BASE);

		//copy
		uint32_t startKey = 0;
		uint32_t endKey = 0xFFFFFFFF;

		uint32_t startLeftPart = (startKey >> BLOCK_BITS);
		uint32_t endLeftPart = (endKey >> BLOCK_BITS);

		//start
		for (uint32_t currHeader = startLeftPart; currHeader <= endLeftPart; currHeader++)
		{
			/*if(currHeader == 941647)
			{
				currHeader = 941647;
			}*/

			HeaderCellInt& headerCell = pHeader[currHeader];
			if (headerCell.Type == 1)
			{
				uint32_t key = (currHeader << BLOCK_BITS) + headerCell.Code;
				uint32_t value = headerCell.Offset;

				pNewArray->insert(key, value);
			}
			else if (headerCell.Type == 2)
			{
				DoublePageInt* pPage = pDoublePages[headerCell.Offset >> 16];

				if (pPage)
				{
					DoubleValueCellInt& valueCell = pPage->pValues[headerCell.Offset & 0xFFFF];

					uint32_t key1 = (currHeader << BLOCK_BITS) + headerCell.Code;
					uint32_t value1 = valueCell.Value1;

					uint32_t key2 = (currHeader << BLOCK_BITS) + valueCell.Code2;
					uint32_t value2 = valueCell.Value2;

					if (headerCell.Code < valueCell.Code2)
					{
						pNewArray->insert(key1, value1);
						pNewArray->insert(key2, value2);
					}
					else
					{
						pNewArray->insert(key2, value2);
						pNewArray->insert(key1, value1);
					}
				}
			}
			else if (headerCell.Type == 3)
			{
				uint32_t startBlockOffset = headerCell.Offset;

				if (currHeader == startLeftPart)
				{
					startBlockOffset += (startKey & BLOCK_SIZE);
				}

				uint32_t endBlockOffset = headerCell.Offset;

				if (currHeader == endLeftPart)
				{
					endBlockOffset += (endKey & BLOCK_SIZE);
				}
				else
				{
					endBlockOffset += BLOCK_SIZE;
				}

				for (uint32_t offset = startBlockOffset, code = 0; offset < endBlockOffset; offset++, code++)
				{
					uint32_t currPage = (offset >> 16);
					MultiplyPageInt* pPage = pMultiplyPages[currPage];

					if (pPage)
					{
						MultiplyValueCellInt& valueCell = pPage->pValues[offset & 0xFFFF];

						if (valueCell.Code == headerCell.Code)
						{
							uint32_t key = (currHeader << BLOCK_BITS) + code;
							uint32_t value = valueCell.Value;

							pNewArray->insert(key, value);
						}
					}
				}
			}
		}

		//destroy
		this->destroy();

		this->pHeader = pNewArray->pHeader;
		this->pDoublePages = pNewArray->pDoublePages;
		this->pMultiplyPages = pNewArray->pMultiplyPages;

		this->CountDoublePage = pNewArray->CountDoublePage;
		this->CountMultiplyPage = pNewArray->CountMultiplyPage;

		this->maxDoubleOffset = pNewArray->maxDoubleOffset;
		this->maxMultiplyOffset = pNewArray->maxMultiplyOffset;

		this->lastHeaderCode = pNewArray->lastHeaderCode;

		delete pNewArray;
	}

	uint32_t getKeysAndValuesByRange(uint32_t* buffer,
		uint32_t size,
		uint32_t startKey,
		uint32_t endKey)
	{
		uint32_t count = 0;

		uint32_t startLeftPart = (startKey >> BLOCK_BITS);
		uint32_t endLeftPart = (endKey >> BLOCK_BITS);

		//start
		for (uint32_t currHeader = startLeftPart; currHeader <= endLeftPart; currHeader++)
		{
			HeaderCellInt& headerCell = pHeader[currHeader];
			if (headerCell.Type == 1)
			{
				uint32_t key = (currHeader << BLOCK_BITS) + headerCell.Code;
				uint32_t value = headerCell.Offset;

				buffer[count++] = key;
				buffer[count++] = value;
			}
			else if (headerCell.Type == 2)
			{
				DoublePageInt* pPage = pDoublePages[headerCell.Offset >> 16];

				if (pPage)
				{
					DoubleValueCellInt& valueCell = pPage->pValues[headerCell.Offset & 0xFFFF];

					uint32_t key1 = (currHeader << BLOCK_BITS) + headerCell.Code;
					uint32_t value1 = valueCell.Value1;

					uint32_t key2 = (currHeader << BLOCK_BITS) + valueCell.Code2;
					uint32_t value2 = valueCell.Value2;

					if (headerCell.Code < valueCell.Code2)
					{
						buffer[count++] = key1;
						buffer[count++] = value1;
						buffer[count++] = key2;
						buffer[count++] = value2;
					}
					else
					{
						buffer[count++] = key2;
						buffer[count++] = value2;
						buffer[count++] = key1;
						buffer[count++] = value1;
					}
				}
			}
			else if (headerCell.Type == 3)
			{
				uint32_t startBlockOffset = headerCell.Offset;

				if (currHeader == startLeftPart)
				{
					startBlockOffset += (startKey & BLOCK_SIZE);
				}

				uint32_t endBlockOffset = headerCell.Offset;

				if (currHeader == endLeftPart)
				{
					endBlockOffset += (endKey & BLOCK_SIZE);
				}
				else
				{
					endBlockOffset += BLOCK_SIZE;
				}

				for (uint32_t offset = startBlockOffset, code = 0; offset <= endBlockOffset; offset++, code++)
				{
					uint32_t currPage = (offset >> 16);
					MultiplyPageInt* pPage = pMultiplyPages[currPage];

					if (pPage)
					{
						MultiplyValueCellInt& valueCell = pPage->pValues[offset & 0xFFFF];

						if (valueCell.Code == headerCell.Code)
						{
							uint32_t key = (currHeader << BLOCK_BITS) + code;
							uint32_t value = valueCell.Value;

							buffer[count++] = key;
							buffer[count++] = value;
						}
					}
				}
			}
		}

		return count;
	}

	uint64_t getTotalMemory()
	{
		return (uint64_t)sizeof(HArrayInt) + HEADER_SIZE * (uint64_t)sizeof(HeaderCellInt) +
			(uint64_t)CountDoublePage * sizeof(DoublePageInt) +
			(uint64_t)CountMultiplyPage * sizeof(MultiplyPageInt);
	}

	void destroy()
	{
		if (pHeader)
		{
			delete[] pHeader;
			pHeader = 0;
		}

		if (pDoublePages)
		{
			for (uint32_t i = 0; i < CountDoublePage; i++)
				delete pDoublePages[i];

			delete[] pDoublePages;
			pDoublePages = 0;
		}

		if (pMultiplyPages)
		{
			for (uint32_t i = 0; i < CountMultiplyPage; i++)
				delete pMultiplyPages[i];

			delete[] pMultiplyPages;
			pMultiplyPages = 0;
		}

		if (pIndexes)
		{
			delete[] pIndexes;
			pIndexes = 0;
		}

		if (pValues)
		{
			delete[] pValues;
			pValues = 0;
		}

		if (pReleaseCells)
		{
			delete[] pReleaseCells;
			pReleaseCells = 0;
		}
	}
};
