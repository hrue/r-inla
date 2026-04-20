# HArray

[![Build status](https://img.shields.io/github/workflow/status/Bazist/HArray/CMake%20Build%20Matrix?style=plastic)](https://img.shields.io/github/workflow/status/Bazist/HArray/CMake%20Build%20Matrix?style=plastic)

## Probably, this is most optimized Trie structure in the World ! Thats all what you need know about this project :)

**HArrayInt** - Key\Value In Memory structure for 32bit keys

**HArray** - Key\Value In Memory structure for keys with variety size

**HArrayChar** - Wrapper for HArray with interfaces: char* key, uint32 keyLen, char* value, uint32 valueLen

**HArrayUniqueIntValueList** - Wrapper for HArray with inteterfaces: uint32* key, uint32 keyLen, List\<uint32\>  listOfUniqueValues

------------------

## [Start overview from Benchmarks](https://github.com/Bazist/HArray/blob/master/Benchmarks.md)

------------------
### Overview

- **High optimized** Trie structure implemented in more than **8K lines**
- Tested on **1 Billion keys** succesfully
- **Keys with variable length** inside one instance of the container
- **Without any Stop World** events such as Rebuild/Rehashing on Insert key.
- **Without any Hash Functions**, the container has adpative algorithm for different nature of keys
- **Scan by Prefix/Scan by Range** functionality as bonus
- **All Keys are sorted**. Ability to set your **Custom Order** of Keys 
- **Predictable** behaviour even in worst case: smoothly consuming memory, almost constantly latency on insert/lookup
- **Prefix Compression** helps reduce memory when keys have the same prefix: urls, file paths, words etc.
- **Serialize/Deserialize** from/to file at a speed close to the max speed of the hard drive
- **Fair Delete** operation with smoothly dismantling each Key. 
  Dead fibres are used for insert new keys, so structure perfect works in massive insert/delete scenarios.

------------------

## Build Library and Benchmarks

### Prerquisites

The following need to be installed and configured on development box:

- C++ development tools.
- [CMake](https://cmake.org/).

### Build on Linux and Mac

```bash
mkdir build
cmake -Bbuild
cd build && make
```

To build release version of the library and benchmarks run instead:

```bash
cmake -Bbuild -DCMAKE_BUILD_TYPE=Release
```

The library (`libharray.a`) and benchmark application (`HArrayBenchmark`) will be created in `build` folder.

### Build on Windows

```cmd
md build
cmake -Bbuild
msbuild build\HArray.sln
```

To build release version run:

```cmd
msbuild build\HArray.sln /property:Configuration=Release
```

The benchmark application `HArrayBenchmark.exe` together with static library `libharray.lib` will be in `build\Debug` or `build\Release` folder.

## Why we love Trie ? Because it has much more functionality and stability than Hashtables and much more faster than Binary Trees. Let's compare properties:

![alt tag](https://raw.githubusercontent.com/Bazist/HArray/master/Images/functionality2.png)

------------------

## Trie ? I heard about Trees and Hastables but don't know anything about Trie
# [Explain me as for Beginners](https://github.com/Bazist/HArray/blob/master/Trie_for_beginners.md)

### Examples

Initialize container

```c++
#include "HArray.h"

HArray ha;
ha.init(); //ha.init(24) set your custom capacity for big datasets
```
Insert a key

```c++
uint32 key[] = { 10, 20, 30, 40 };
uint32 keyLen = sizeof(key) / 4; //in key segments
uint32 value = 1;

ha.insert(key, keyLen, value);
```

Get value by key. Will return 0 if key is not found

```c++
uint32 value;
if(ha.getValueByKey(key, keyLen, value))
{
   printf("%d", value)
}
else
{
   //key is not found
}
```

Get all keys by range from min key to max key. 

```c++
HArrayPair pairs[5];
uint32 pairsSize = 5;

uint32 minKey[] = { 10, 10 };
uint32 minKeyLen = sizeof(minKey) / 4; //in key segments

uint32 maxKey[] = { 20, 20 };
uint32 maxKeyLen = sizeof(maxKey) / 4; //in key segments

uint32 count = ha.getKeysAndValuesByRange(pairs, pairsSize, minKey, minKeyLen, maxKey, maxKeyLen);
for (uint32 i = 0; i < count; i++)
{
   uint32* key = pairs[i].Key;
   uint32 keyLen = pairs[i].KeyLen;

   uint32 value = pairs[i].Value;
   
   //here your code
}
```

Scan all keys through visitor

```c++
bool visitor(uint32* key, uint32 keyLen, uint32 value, uchar8 valueType, void* pData)
{
   //handle founded key here
   // ...

   //return true to continue scaning or false otherwise
   return true;
}

//Start scanning

void* pData = 0; //useful data in context

ha.scanKeysAndValues(&visitor, pData);
```

Serialize/Deserialize containter from/to file

```c++
ha.saveToFile("c:\\dump");

ha.loadFromFile("c:\\dump");
```

Check if container has part of key

```c++
uint32 key[] = { 10, 20, 30 };
uint32 keyLen = sizeof(key) / 4; //in key segments

if(ha.hasPartKey(key, keyLen))
{
   //code here
}
```

Set specific comparator to redefine order of keys in the container.

```c++
float key[] = { 10.0, 20.0, 30.0 };
uint32 keyLen = sizeof(key) / 4; //in key segments

uint32 value = 1;

//Set float comparator for right sorting
//Another options: setStrComparator, setInt32Comparator, setUInt32Comparator 
//or define your custom comparator through setCustomComparator
ha.setFloatComparator();

ha.insert((uint32*)key, keyLen, value);

```

Delete Key and Value from container

```c++
ha.delValueByKey(key, keyLen);
```

### How to Run

```c++
git clone github.com/Bazist/HArray.git
cd HArray/HArray
make
./HArray
```

### License

The code is licensed under the GNU General Public License v3.0, see the [LICENSE file](LICENSE) for details.
