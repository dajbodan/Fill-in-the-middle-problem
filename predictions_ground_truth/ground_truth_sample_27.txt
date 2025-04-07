#ifndef __PROGTEST__
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <climits>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include <set>
#include <queue>
#include <memory>
#include <functional>
#include <stdexcept>
using namespace std;
#endif /* __PROGTEST__ */

enum GRAPH
{
  LEAF = 0, INNODE = 1, DEFAULT = 2
};


typedef struct tUTF8
{
  enum Bytes
  {
      UNKNOWN, ONE, TWO, THREE, FOUR
  };

    unsigned char data[4];
    Bytes lengthData;
    tUTF8(void) 
    {
      for(int i = 0; i < 4; i++)
        data[i] = '\0'; 
      lengthData = UNKNOWN;
      }
    tUTF8(char * inputData, Bytes inputLengthData) : lengthData(inputLengthData) 
    {
      for(int i = 0; i < inputLengthData; i++)
        data[i] = inputData[i];
    }
    bool operator < (const tUTF8 & rhs) const
    {
        uint64_t dataInIntLHS = static_cast<unsigned char>(this->data[0]) * 1 
                                + static_cast<unsigned char>(this->data[1]) * 256 
                                + static_cast<unsigned char>(this->data[2]) * 65536
                                + static_cast<unsigned char>(this->data[3]) * 16777216;
        uint64_t dataInIntRHS = static_cast<unsigned char>(rhs.data[0]) * 1
                                + static_cast<unsigned char>(rhs.data[1]) * 256
                                + static_cast<unsigned char>(rhs.data[2]) * 65536
                                + static_cast<unsigned char>(rhs.data[3]) * 16777216;
        return dataInIntLHS < dataInIntRHS;
    } 
}tUTF8;

typedef struct tNodeCompression
{
  tUTF8 utf8;
  int frequency;
  GRAPH typeNode;
  tNodeCompression()
  {
    frequency = 0;
    typeNode = GRAPH :: DEFAULT;
  }
  tNodeCompression(tUTF8 inputUTF8) 
  {
    frequency = 1;
    utf8 = inputUTF8;
    typeNode = GRAPH :: DEFAULT;
  }
  tNodeCompression(tUTF8 inputUTF8, int inputFreq, GRAPH inputTypeNode )
  {
    frequency = inputFreq;
    utf8 = inputUTF8;
    typeNode = inputTypeNode;
  }
  
}tNodeCompression;


typedef struct BinTreeCompression
{
    tNodeCompression  data;
    BinTreeCompression * left;
    BinTreeCompression * right;

    BinTreeCompression()
    {
      left = nullptr;
      right = nullptr;
    }
    BinTreeCompression(GRAPH inputTypeNode) : left(nullptr), right(nullptr)  {}
    BinTreeCompression(
      BinTreeCompression * inputLeft,
      BinTreeCompression * inputRight,
      tNodeCompression inputData,
      GRAPH inputTypeNode)
    {
      left = inputLeft;
      right = inputRight;
      for(int i = 0; i < inputData.utf8.lengthData; i++)
        data.utf8.data[i]   = inputData.utf8.data[i];
      data.frequency        = inputData.frequency;
      data.utf8.lengthData  = inputData.utf8.lengthData;
      data.typeNode         = inputTypeNode;
    }
    BinTreeCompression(BinTreeCompression * inputLeft, BinTreeCompression * inputRight, int frequency)
    {
      data.typeNode = GRAPH ::INNODE;
      data.frequency = frequency;
      left = inputLeft;
      right = inputRight;
    }
    friend bool operator == (const BinTreeCompression * lhs, const tNodeCompression & rhs);
    
    
}BinTreeCompression;

bool operator == (const BinTreeCompression * lhs, const tNodeCompression & rhs)
{
  if(lhs->data.frequency == rhs.frequency)
    return true;
  return false;
}


enum STATES
{
    readingINNODE, readingLIST, readingUTF8, UNDEFINED
};

typedef struct BinTreeNode
{
  BinTreeNode * left;
  BinTreeNode * right;

  GRAPH typeNode;

  tUTF8 utf8;

  BinTreeNode() : left(nullptr), right(nullptr), typeNode(GRAPH :: DEFAULT){}
}BinTreeNode;


typedef struct tReadOneByte
{
  char oneByte[1];
  int readIndex;
  tReadOneByte() : readIndex(0) { oneByte[0] = '\0'; }
}tReadOneByte;

class myCompratatorForCompression
{
  public : 
    int operator() (const tNodeCompression & p1, const tNodeCompression & p2)
    {
      return p1.frequency > p2.frequency; 
    }
};

class HuffmanTree
{
  private : 
    BinTreeNode * root;
    int cntPossibleLeafs;

  public : 
    HuffmanTree(void) : cntPossibleLeafs(2) { root = nullptr; } // leafs actualy is 0, but after creation of root, it can be != 0
    HuffmanTree(GRAPH typeNode);
    
    void insertValue( GRAPH typeNode, tUTF8 utf8 = tUTF8());
    int getCntPossibleLeafs() { return cntPossibleLeafs; }
    BinTreeNode * getRoot() { return root; }
    BinTreeNode * getLeftNode(BinTreeNode * node) { if(node != nullptr) return node->left; return nullptr;}
    BinTreeNode * getRightNode(BinTreeNode * node) { if(node != nullptr) return node->right; return nullptr;}



    ~HuffmanTree(void);
    private :
      static void recDeallocateTree(BinTreeNode * node);
      BinTreeNode * recFindVal(BinTreeNode * node);
      BinTreeNode * createVal(GRAPH typeNode, tUTF8 utf8 = tUTF8());

};

class HuffmanCodeDecode
{
  HuffmanTree * tree;

  public : 

    HuffmanCodeDecode(void);
    ~HuffmanCodeDecode(void) { delete tree; }
    HuffmanTree * getTree(void) { return tree; }
};

HuffmanCodeDecode :: HuffmanCodeDecode(void) { tree = new HuffmanTree; }
HuffmanTree :: HuffmanTree(GRAPH typeNode) : cntPossibleLeafs(0)
{
  root = createVal(typeNode);
  typeNode = GRAPH :: DEFAULT;
}

HuffmanTree :: ~HuffmanTree(void) { recDeallocateTree(root); }

BinTreeNode * HuffmanTree :: recFindVal(BinTreeNode * node)
{
  BinTreeNode * res = node;
  if(node->typeNode == GRAPH ::LEAF)
    return nullptr;
  //if(node ->typeNode ) // NEED TO THINK, IF ROOT CAN BE LEAF
  if(node -> left == nullptr)
    return node;

  if(node -> left != nullptr)
    res = recFindVal(node -> left);
  if(res != nullptr)
      return res;
  if(node -> right == nullptr)
    return node;

  if(node -> right != nullptr)
    res = recFindVal(node -> right);
  if(res != nullptr)
    return res;
  else 
    return nullptr;
  
}


void HuffmanTree :: recDeallocateTree(BinTreeNode * node)
{
  if(node == nullptr)
    return ;
  if(node -> left != nullptr)
    recDeallocateTree(node -> left);
  if(node -> right != nullptr)
    recDeallocateTree(node -> right);
  delete node;
}

void HuffmanTree :: insertValue( GRAPH typeNode, tUTF8 utf8)
{
  if(root == nullptr)
  {
    if(typeNode == GRAPH :: LEAF)
      createVal(typeNode, utf8);
    else if(typeNode == GRAPH :: INNODE)
      root = createVal(typeNode);
    cntPossibleLeafs = 2;
    return ;
  }
  BinTreeNode * nodeInsert = recFindVal(root);
  if(nodeInsert -> left == nullptr || nodeInsert->right == nullptr)
    cntPossibleLeafs--;

  if(typeNode == GRAPH :: LEAF)
  {
    if(nodeInsert -> left == nullptr)
      nodeInsert -> left = createVal(typeNode, utf8);
    else if( nodeInsert -> right == nullptr)
      nodeInsert -> right = createVal(typeNode, utf8);
  }
  else if(typeNode == GRAPH :: INNODE)
  {
    if(nodeInsert -> left == nullptr)
      nodeInsert -> left = createVal(typeNode);
    else if( nodeInsert -> right == nullptr)
      nodeInsert -> right = createVal(typeNode);
  }
}

BinTreeNode * HuffmanTree :: createVal(GRAPH typeNode, tUTF8 utf8 )
{
  BinTreeNode * node = new BinTreeNode;
  node -> left = nullptr;
  node -> right = nullptr;
  node -> typeNode = typeNode;
  if(typeNode == GRAPH :: LEAF)
  {
    node -> utf8 = utf8;
  }
  else if(typeNode == GRAPH ::INNODE) 
  {
    cntPossibleLeafs += 2;
  }
  return node;
}

unsigned char changeReadDirection(int & i, std :: queue <tReadOneByte> & readByteQueue)
{
  i = 0;
  readByteQueue.pop();
  return readByteQueue.front().oneByte[0];
}

HuffmanTree * readBitsInHeader
    ( 
      HuffmanTree * tree,
      std :: queue <tReadOneByte> & readByteQueue,
      STATES & state,
      const int lengthByte = 8
    )
{

    unsigned char byte = readByteQueue.front().oneByte[0];
    const char MSB = 1 << (lengthByte - 1);
    byte <<= readByteQueue.front().readIndex;
    for(int i = readByteQueue.front().readIndex ; state != readingUTF8; i++)
    {
        if(i == 8)
          byte = changeReadDirection(i,readByteQueue);
        

        state = ((byte & MSB) == 0b0) ? readingINNODE : readingLIST ; // NEXT STATE
        switch(state)
        {
          case readingINNODE :
            
            tree -> insertValue(GRAPH :: INNODE);
            break;
          case readingLIST : 
           // tree -> insertValue(GRAPH :: LEAF);
            readByteQueue.front().readIndex = i + 1; // next Index
            if(readByteQueue.front().readIndex == 8)
              byte = changeReadDirection(i, readByteQueue);
            state = readingUTF8;
            return tree;
          default :
            break;
        }
        byte <<= 1;

    }
  return tree;
}
const int lengthFour = 4;

void makeByte(std :: queue<tReadOneByte> & readByteQueue, unsigned char data[lengthFour], int where, const int lengthBytes = 8)
{
  unsigned char res = readByteQueue.front().oneByte[0] << readByteQueue.front().readIndex;
  int tempPrevLength = readByteQueue.front().readIndex; 
  readByteQueue.pop();
  unsigned char temp = readByteQueue.front().oneByte[0];
  temp >>= lengthBytes - tempPrevLength;
  res |= temp;
  readByteQueue.front().readIndex = tempPrevLength;
  data[where] = res;
}

void copyBytes(char * d, unsigned char * s, int length)
{
  for(int i = 0; i < length; i++)
    d[i] = s[i];
}

bool isLegitUTF8(unsigned char s[4], const int lengthByteInBit = 8)
{
  unsigned char tempS[4];
  for(int i = 0; i < 4; i++)
    tempS[i] = s[i];
  unsigned char maxUTF8[4] = {0xF4, 0x8F, 0xBF, 0xBF };
  for(int i = 0; i < 4; i ++)
  {
    for(int j = 0; j < 8; j++, maxUTF8[i] <<= 1, tempS[i] <<= 1)
    {
        if(((maxUTF8[i] & 0x80) ^ (tempS[i] & 0x80)) != 0b0)
          return ((maxUTF8[i] & 0x80) != 0b0
                  && (tempS[i] & 0x80) == 0b0); 
          
    }
  }
  return true;
}



void readUTF8
        (
          std :: queue<tReadOneByte> & readByteQueue,
          tUTF8 & utf8,
          const int lengthByte = 8 
        )
{
  unsigned char readFourBytes[lengthFour];
  makeByte(readByteQueue,readFourBytes, 0);
  if((readFourBytes[0] & 0x80) == 0x00 )
  {
    copyBytes((char *)utf8.data, readFourBytes, 1);
    utf8.lengthData = tUTF8 :: Bytes :: ONE;
    return ;
  }
  
  makeByte(readByteQueue, readFourBytes, 1);
  if((readFourBytes[0] & 0xE0) == 0xC0 && (readFourBytes[1] & 0xC0) == 0x80)
  {
    copyBytes((char *)utf8.data, readFourBytes, 2);
    utf8.lengthData = tUTF8 :: Bytes :: TWO;
    return ;
  }
   
  
  makeByte(readByteQueue, readFourBytes, 2);
  if(((readFourBytes[0] & 0xF0) == 0xE0) && ((readFourBytes[1] & 0xC0) == 0x80) && ((readFourBytes[2] & 0xC0) == 0x80)) 
  {
    copyBytes((char *)utf8.data, readFourBytes, 3);
    utf8.lengthData = tUTF8 :: Bytes :: THREE;
    return ;
  }
  
  
  makeByte(readByteQueue, readFourBytes, 3);
  if(((readFourBytes[0] & 0xF8) == 0xF0) 
      && ((readFourBytes[1] & 0xC0) == 0x80)
      && ((readFourBytes[2] & 0xC0) == 0x80)
      && ((readFourBytes[3] & 0xC0) == 0x80)
      && isLegitUTF8(readFourBytes)
    )
    {
      copyBytes((char *)utf8.data, readFourBytes, 4);
      utf8.lengthData = tUTF8 :: Bytes :: FOUR;
    return ;
    }
    
  
  //std :: cout << "size : " <<  readByteQueue.size() << "\n";
  throw std :: exception();
}

void readHeader
        (
          HuffmanTree * tree,
          std :: queue<tReadOneByte> & readByteQueue  
        )
{

  STATES state = STATES :: readingINNODE;
   while(tree->getCntPossibleLeafs() != 0)
  {
      
      tree = readBitsInHeader(tree, readByteQueue ,state);
      if(state == STATES :: readingUTF8)
      {
        tUTF8 utf8;
        readUTF8( readByteQueue, utf8);
        tree->insertValue(GRAPH :: LEAF, utf8);
        state = STATES :: UNDEFINED;
      }
       
  }
}



bool readInputFile(std :: queue<tReadOneByte> & readQueue, const char * fileName)
{
  std :: ifstream inputFile;
  
  inputFile.open(fileName, std :: ios_base :: binary);
  if(! inputFile.is_open())
    return true;
  tReadOneByte varReadByte;
  const int byte = 1;
  while(true)
  {
    inputFile.read(varReadByte.oneByte,byte);
    if(inputFile.eof())
      break;
    varReadByte.readIndex = 0;
    readQueue.push(varReadByte);
  }
  if(inputFile.bad())
    return true;
  inputFile.close();
  return false;
}

bool isLastChunk(char source, int & startRead, const int lengthByte = 8)
{
  unsigned char temp = source << startRead;
  unsigned char mask = 0x80;
  startRead++;
  return (temp & mask ) == (unsigned char)0b0;
}


int readSize(
        std :: queue<tReadOneByte> & readQueue,
        const int lengthSizeInBits = 12) // need to fix. What if start read Index would be at 7
{
  int res = 0;
  const int sizeByteInBits = 8;
  tReadOneByte arrReadSize[2];
  arrReadSize[0] = readQueue.front();
  readQueue.pop();
  arrReadSize[0].oneByte[0] <<= arrReadSize[0].readIndex;

  int mod = 1 << (lengthSizeInBits - 1);
  int maskMostLeftBit = 0x80;
  int readBits = 0;

  for(int i = 0; i < lengthSizeInBits; i++, arrReadSize[0].oneByte[0] <<= 1, mod >>= 1, readBits++)
  {
    if(arrReadSize[0].readIndex + readBits == sizeByteInBits)
    {
      if(readBits == sizeByteInBits)
        readQueue.pop();
      arrReadSize[1] = readQueue.front();
      arrReadSize[0] = arrReadSize[1];
      readBits = 0;
      
    }
    if((arrReadSize[0].oneByte[0] & maskMostLeftBit) != 0b0)
      res += mod;
  }
  readQueue.front().readIndex = readBits;
  return res;

}

tUTF8 recReadTree
(
  
  BinTreeNode * node,
  std :: queue<tReadOneByte> & readQueue,
  const int & sizeByteIntBits = 8
)
{
  if(readQueue.front().readIndex == sizeByteIntBits)
  {
    readQueue.pop();
  }
  
  if(node->typeNode == GRAPH :: LEAF)
    return node->utf8;
  
  tReadOneByte readVal = readQueue.front();//readQueue.front().oneByte[0] << readQueue.front().readIndex;
  readVal.oneByte[0] <<= readVal.readIndex;

  if((readVal.oneByte[0] & (unsigned char)0x80) == 0b0)
  {
    readQueue.front().readIndex++;
    return recReadTree( node->left, readQueue);
  }
  else
  {
    readQueue.front().readIndex++;
    return recReadTree(node->right, readQueue);
  }
}

void copy(char d[4], tUTF8 & s)
{
  for(int i = 0; i < s.lengthData; i++)
    d[i] = s.data[i];
}

void readChunk(std :: queue<tUTF8> & out, BinTreeNode * root, std :: queue<tReadOneByte> & readQueue, int size)
{
    for(int i = 0; i < size; i++)
    {
      if(readQueue.size() == 0)
      {
          throw std :: exception();
      }
      tUTF8 res = recReadTree(root, readQueue);
      
      
      char tempRes[4];
      copy(tempRes, res);
      out.push(res);
    }
}

void readLastChunk
(
  std :: queue<tUTF8> & out,
  BinTreeNode * root,
  std :: queue<tReadOneByte> & readQueue
)
{
  int size = readSize(readQueue);
  readChunk(out, root, readQueue, size);

}

bool readCommpresedData
(
  std :: queue<tUTF8> & out,
  BinTreeNode * root,
  std :: queue<tReadOneByte> & readQueue
  
)
{
  
  tReadOneByte * actualReadElement = &readQueue.front();
  //int some = 0;
  while(! isLastChunk(actualReadElement->oneByte[0], actualReadElement->readIndex))
  {
    const int lengthChunk = 4096;
    readChunk(out, root, readQueue,lengthChunk );
    if(readQueue.front().readIndex == 8)
    {
      readQueue.pop();
    }
    actualReadElement = & readQueue.front();
  }
  readLastChunk(out, root, readQueue);
  return false;
}



bool writeOutput(const char * outFileName, std :: queue<tUTF8> & outQueue)
{
  std :: ofstream out;
  out.open(outFileName, std :: ios_base :: binary);
  if(!out.is_open() || !out)
  {
   // out.close();
    return true;
  }
  while(outQueue.size() != 0)
  {
    tUTF8 res = outQueue.front();
    outQueue.pop();
    char tempRes[4];
    copy(tempRes, res);
    out.write(tempRes,res.lengthData);
  }
  if(out.bad())
    return true;
  out.close();
  return false;
}

bool decompressFile ( const char * inFileName, const char * outFileName )
{
  HuffmanCodeDecode huffmanCode;
  std :: queue<tReadOneByte> readQueue;
  if(readInputFile(readQueue, inFileName))
  {
    return false;
  }
  std :: queue<tUTF8> out; 
  try
  {
    readHeader(huffmanCode.getTree(), readQueue);
  }
  catch(std :: exception & ex)
  {
    return false;
  }
  try{
    if(readCommpresedData(out, huffmanCode.getTree()->getRoot(), readQueue ))
      return false;
  }
   catch(std :: exception & ex)
  {
    return false;
  }
  if(writeOutput(outFileName,out))
    return false;
  return true;
}

void printCharacter(std :: vector<tNodeCompression> & vecCompres)
{
  for(tNodeCompression & temp : vecCompres)
  {
    std :: cout << "Char : ";
    for(int i = 0; i < temp.utf8.lengthData; i++)
        std :: cout << temp.utf8.data[i];
    std :: cout << " , Freq : " << temp.frequency << "\n";  
  }
}

void printCharacter(tNodeCompression & s)
{
    std :: cout << "Char : ";
    for(int i = 0; i < s.utf8.lengthData; i++)
        std :: cout << s.utf8.data[i];
    std :: cout << " , Freq : " << s.frequency << "\n"; 
}

bool isSameUTF8(tUTF8 & first, const tUTF8 second)
{
  if(first.lengthData != second.lengthData)
    return false;
  for(int i = 0; i < first.lengthData; i++)
      if(first.data[i] != second.data[i])
        return false;
  return true;
}

bool isPresent(std :: vector<tNodeCompression> & source, const tUTF8 what)
{
  for(tNodeCompression & temp : source)
  {
    if( isSameUTF8(temp.utf8, what))
    {
      temp.frequency++;
      return true;
    }
  }
  return false;
}

void calcFrequency(const std :: vector<tUTF8> & vecOnlyLetters, std :: vector<tNodeCompression> & vecFrequency )
{
    for(const tUTF8 & letter : vecOnlyLetters)
    {
        if(! isPresent(vecFrequency, letter))
        {
            vecFrequency.push_back(tNodeCompression(letter, 1, GRAPH :: LEAF));
        }
    }
}

void makePQForCommpresion(
  std :: priority_queue<tNodeCompression, std :: vector<tNodeCompression>, myCompratatorForCompression> & pq,
  std :: vector<tNodeCompression> & vecFrequency
  )
{
  for(size_t i = 0; i < vecFrequency.size(); i++)
    pq.push(vecFrequency.at(i));
}

BinTreeCompression * findNode(std :: vector<BinTreeCompression * > & vecPossibleNodes, tNodeCompression & cmp)
{
  
  for(size_t i = 0; i < vecPossibleNodes.size(); i++ )
  {
    if(vecPossibleNodes.at(i) == cmp)
    {
      BinTreeCompression * temp = vecPossibleNodes.at(i);
      vecPossibleNodes.erase(vecPossibleNodes.begin() + i);
      return temp;
    }
  }
  throw  std :: exception(); 
}

BinTreeCompression * makeAHuffmanTree( std :: vector<tNodeCompression> & vecFrequency)
{
  BinTreeCompression * root = nullptr;
  
  std :: priority_queue<tNodeCompression, std :: vector<tNodeCompression>, myCompratatorForCompression > pq;
  std :: vector<BinTreeCompression * >  vecPossibleNodes;
  makePQForCommpresion(pq, vecFrequency);
  while(pq.size() != 1)
  {
    tNodeCompression leftNode = pq.top();
    pq.pop();
    tNodeCompression rightNode = pq.top();
    pq.pop();

    if(leftNode.typeNode == GRAPH :: LEAF)
    {
      vecPossibleNodes.push_back(new BinTreeCompression(nullptr,nullptr, leftNode, GRAPH :: LEAF)); 
    }
    if(rightNode.typeNode == GRAPH :: LEAF)
    {
      vecPossibleNodes.push_back(new BinTreeCompression(nullptr, nullptr, rightNode, GRAPH :: LEAF));
    }
    BinTreeCompression * left = findNode(vecPossibleNodes, leftNode);
    BinTreeCompression * right = findNode(vecPossibleNodes, rightNode);
    BinTreeCompression * insertNode = new BinTreeCompression(left, right, left->data.frequency + right->data.frequency);
    vecPossibleNodes.push_back(insertNode);
    pq.push(insertNode->data);
  }
  root = vecPossibleNodes.at(0);
  return root;
}

void readLetterFromInputFile(std :: vector<tUTF8> & vecOnlyLetters, std :: queue<tReadOneByte> & readQueue)
{
   while(readQueue.size() != 0)
    {
      tUTF8 temp;
      readUTF8(readQueue,temp);
      vecOnlyLetters.push_back(temp);
    }
}

void writeLetterToOutput(std :: deque<tReadOneByte> & outputQueue, const tUTF8 & utf8)
{
  
  for(int i = 0; i < utf8.lengthData; i++)
  {
    char temp = utf8.data[i];
    for(int j = 0; j < 8; j++, outputQueue.back().readIndex++, temp <<= 1)
    {
        if(outputQueue.back().readIndex == 8)
        {
          outputQueue.push_back(tReadOneByte());
        }
        
        outputQueue.back().oneByte[0] |= ((temp & 0x80) >> outputQueue.back().readIndex);
    }

  }
}

void recOutputTree(std :: deque<tReadOneByte> & outputQueue, BinTreeCompression * root, const int lengthBytesInBits = 8)
{
  if(outputQueue.back().readIndex == 8) // was 8, was working 
  {
    outputQueue.push_back(tReadOneByte());
  } 
  if(root->data.typeNode == GRAPH :: INNODE)
  {
    outputQueue.back().readIndex++;
    outputQueue.back().oneByte[0] |= 0 << ( lengthBytesInBits -  outputQueue.back().readIndex);
    
  }
  if(root->data.typeNode == GRAPH :: LEAF)
  {
    outputQueue.back().readIndex++;
    outputQueue.back().oneByte[0] |= 1 << (lengthBytesInBits - outputQueue.back().readIndex);
    writeLetterToOutput(outputQueue, root->data.utf8);
    
  }
  if(root->left != nullptr)
    recOutputTree(outputQueue, root->left);

  if(root->right != nullptr)
    recOutputTree(outputQueue, root->right);
}

bool writeOutput(const char * outFileName, std :: deque<tReadOneByte> & outputDeq)
{
  std :: ofstream out;
  out.open(outFileName, std :: ios_base :: binary);
  if(!out.is_open() || !out)
  {
   // out.close();
    return true;
  }
  while(outputDeq.size() != 0)
  {
    char res[2];
    res[0] = outputDeq.front().oneByte[0];
    outputDeq.pop_front();
    out.write(res,1);
  }
  if(out.bad())
    return true;
  out.close();
  return false; 
}

bool isSame(tUTF8 &firstCMP, tUTF8 &secondCMP)
{
  if(firstCMP.lengthData != secondCMP.lengthData)
    return false;
  for(int i = 0; i < firstCMP.lengthData; i++)
    if(firstCMP.data[i] != secondCMP.data[i])
      return false;
  return true;
}

void makeMask(unsigned char & mask, size_t readIndex, const int lengthByteInBits = 8)
{
  for(size_t i = 0; i < readIndex; i++)
  {
    mask |= 1 << (lengthByteInBits - i);
  }
}

void checkIndex(std :: vector<tReadOneByte> & res)
{
  if(res.back().readIndex == 8)
  {
    res.push_back(tReadOneByte());
  }
}

void deallocateCompressedTree(BinTreeCompression * root)
{
  if(root->left != nullptr)
    deallocateCompressedTree(root->left);
  if(root->right != nullptr)
    deallocateCompressedTree(root->right);
  delete root;
}

std :: vector<tReadOneByte> & recTranslateLetterFromUnicodeToHuffman
(
  tUTF8 & source,
  std :: vector<tReadOneByte> & res,
  BinTreeCompression * root,
  bool & found,
  const int lengthByteInBits = 8
)
{
  checkIndex(res);
  if(root->data.typeNode == GRAPH :: LEAF)
  {
    if(isSame(source, root->data.utf8))
    {
        found = true;
        return res;
    }
  }
  if(root->left != nullptr)
  {
    res.back().readIndex++;
    
    unsigned char mask = '\0';
    makeMask(mask, res.back().readIndex);
    res.back().oneByte[0] &= mask; //0 << (lengthByteInBits - res.at(pos).readIndex);
    res = recTranslateLetterFromUnicodeToHuffman(source,res, root->left, found);
    if(found == true)
      return res;
    res.back().readIndex--;
    
  }
  if(root->right != nullptr)
  {
    
    res.back().readIndex++;
    
    unsigned char mask = '\0';
    makeMask(mask, res.back().readIndex);
    res.back().oneByte[0] &= mask;
    res.back().oneByte[0] |= 1 << (lengthByteInBits - res.back().readIndex);
    res = recTranslateLetterFromUnicodeToHuffman(source, res, root->right, found);
    if(found == true)
      return res;
    res.back().readIndex--;
  }
  if( res.size() != 1 && res.back().readIndex == 0)
  {
    res.pop_back();
    //res.at(pos).readIndex--;
  }
  return res;
}


void translateLengthToBits(std :: deque<tReadOneByte> & outputDeq, int size, const int lengthByteInBits = 8)
{
  int mod = 2048;
  for(int i = 0; i < 12; i++, mod /= 2)
  {
    if(outputDeq.back().readIndex == 8)
      outputDeq.push_back(tReadOneByte());
    outputDeq.back().readIndex++;
    if(size / mod != 0)
    {
      outputDeq.back().oneByte[0] |= (1 << (lengthByteInBits - outputDeq.back().readIndex) );
      size -= mod;
    }
    else 
      outputDeq.back().oneByte[0] |= (0 << (lengthByteInBits - outputDeq.back().readIndex));
  }
}

void makeHuffmanCode
          (
            std :: deque<tReadOneByte> & outputDeq,
            std :: vector<tUTF8> & vecOnlyLetters,
            std :: map<tUTF8, std :: vector<tReadOneByte>> & mapLetterIntoHuffmanCode
          )
{
  std :: vector<tReadOneByte> letterInHuffmanCode = mapLetterIntoHuffmanCode[vecOnlyLetters.front()]; //recTranslateLetterFromUnicodeToHuffman(vecOnlyLetters.front(), letterInHuffmanCode, root,  found );
  for(size_t i = 0; i < letterInHuffmanCode.size(); i++)
    {
      unsigned char temp = static_cast<unsigned char>(letterInHuffmanCode.at(i).oneByte[0]);
      for(int j = 0; j < letterInHuffmanCode.at(i).readIndex; j++, outputDeq.back().readIndex++ ,temp <<= 1)
      {
        if(outputDeq.back().readIndex == 8)
          outputDeq.push_back(tReadOneByte());
        
        outputDeq.back().oneByte[0] |= ((temp & 0x80) >> outputDeq.back().readIndex);
      }
    }
    vecOnlyLetters.erase(vecOnlyLetters.begin());
} 

void makeChunks(
        std :: vector<tUTF8> & vecOnlyLetters,
        std :: deque<tReadOneByte> & outputDeq,
        std :: map<tUTF8, std :: vector<tReadOneByte>> & mapLetterIntoHuffmanCode,
        const int lengthByteInBits = 8)
{
  

  while(vecOnlyLetters.size() > 4095)
  {
    //std :: cout << "Compressed :: makeChunks :: localVar :: outputDeq, Iterate : " << outputDeq.size() << "\n";
    if(outputDeq.back().readIndex == 8)
      outputDeq.push_back(tReadOneByte());
    outputDeq.back().readIndex++;
   
    //std :: cout << "Iterate: Val :" << outputDeq.back().readIndex << "\n";
    outputDeq.back().oneByte[0] |= 1 << (lengthByteInBits - outputDeq.back().readIndex);
    for(int i = 0; i < 4096; i++)
      makeHuffmanCode(outputDeq,vecOnlyLetters, mapLetterIntoHuffmanCode  );
  }
  outputDeq.back().readIndex++;
  outputDeq.back().oneByte[0] |= 0 << (lengthByteInBits - outputDeq.back().readIndex);
  //std :: cout << "Compressed :: makeChunks :: localVar :: outputDeq " << outputDeq.size() << "\n";
  //std :: cout << "Compressed :: makeChunks :: localVar :: vecOnlyLetters " << vecOnlyLetters.size() << "\n";
  translateLengthToBits(outputDeq,vecOnlyLetters.size());
  size_t size = vecOnlyLetters.size();
  for(size_t i = 0; i < size; i++)
  {
    makeHuffmanCode(outputDeq,vecOnlyLetters, mapLetterIntoHuffmanCode);
  }
   // std :: cout << "Compressed :: makeChunks :: localVar :: outputDeq " << outputDeq.size() << "\n";

}


void recMakeMap
          (
            std :: map<tUTF8, std :: vector<tReadOneByte>> & map,
            BinTreeCompression * root,
            std :: vector<tReadOneByte> & insertVal,
            const int lengthByteInBits = 8
          )
{
  if(root->data.typeNode == GRAPH :: LEAF)
  {
    map.insert( std :: pair<tUTF8, std :: vector<tReadOneByte>>(root->data.utf8, insertVal));
    return ;
  }
  checkIndex(insertVal);
  if(root->left != nullptr)
  {
    insertVal.back().readIndex++;
    unsigned char mask = '\0';
    makeMask(mask, insertVal.back().readIndex);
    insertVal.back().oneByte[0] &= mask;
    recMakeMap(map, root->left, insertVal);
    insertVal.back().readIndex--;
  }
  if(root->right != nullptr)
  {
    insertVal.back().readIndex++;
    unsigned char mask = '\0';
    makeMask(mask, insertVal.back().readIndex);
    insertVal.back().oneByte[0] &= mask;
    insertVal.back().oneByte[0] |= 1 << (lengthByteInBits - insertVal.back().readIndex);
    recMakeMap(map, root->right, insertVal);
    insertVal.back().readIndex--;
  }
  if(insertVal.size() != 1 && insertVal.back().readIndex == 0)
  {
    insertVal.pop_back();
  }
}

bool compressFile ( const char * inFileName, const char * outFileName )
{
  std :: queue<tReadOneByte> readQueue;
  std :: vector<tUTF8> vecOnlyLetters;
  if(readInputFile(readQueue,inFileName))
    return false;
  try 
  {
    readLetterFromInputFile(vecOnlyLetters, readQueue);
  }
  catch(std :: exception & ex)
  {
    return false;
  }

  std :: vector<tNodeCompression> vecFrequency;
  
  calcFrequency(vecOnlyLetters,vecFrequency);

  //printCharacter(vecFrequency);
  BinTreeCompression * root =  makeAHuffmanTree(vecFrequency); 

  std :: deque<tReadOneByte> outputDeq;
  outputDeq.push_back(tReadOneByte()); // initializeVal
  std :: map<tUTF8, std :: vector<tReadOneByte>> mapLetterIntoHuffmanCode; 
  recOutputTree(outputDeq, root);
  std :: vector<tReadOneByte> tempInsertVal;
  tempInsertVal.push_back(tReadOneByte());
  recMakeMap(mapLetterIntoHuffmanCode , root, tempInsertVal);
  makeChunks(vecOnlyLetters,outputDeq, mapLetterIntoHuffmanCode );
  if(writeOutput(outFileName,outputDeq))
  {
    deallocateCompressedTree(root);
    return false;
  }
  deallocateCompressedTree(root);
  return true;
}
#ifndef __PROGTEST__
bool identicalFiles ( const char * fileName1, const char * fileName2 )
{
  std :: ifstream firstFile;
  std :: ifstream secondFile;
  firstFile.open(fileName1, std :: ios_base :: binary);
  secondFile.open(fileName2, std :: ios_base :: binary);
  const int bytes4 = 4;
  char first[4];
  char second[4];
  while( true )
  {
    firstFile.read(first, bytes4);
    secondFile.read(second, bytes4);
    if(firstFile.eof() || secondFile.eof())
        break;  
    for(int i = 0; i < bytes4; i++)
      if(first[i] != second[i])
      {
        //std :: cout << "Where : " << firstFile.tellg() << "\n";
        return false;
      }
  }
  
  return  firstFile.eof() && secondFile.eof(); 
}

int main ( void )
{

  /*assert(compressFile("some.txt", "tempfile"));
  assert(decompressFile("tempfile", "tempfile2"));
  assert( identicalFiles("tempfile2", "some.txt"));*/

  /*assert(compressFile("tests/test0.orig", "tempfile"));
  assert(decompressFile("tempfile", "tempfile2"));
  assert( identicalFiles("tempfile2", "tests/test0.orig"));
  
  assert ( compressFile ( "tests/test1.orig", "tempfile" ) );
  assert ( decompressFile ( "tempFile", "tempfile2" ) );
  assert( identicalFiles("tempfile2", "tests/test1.orig"));

  assert ( compressFile ( "tests/test2.orig", "tempfile" ) );
  assert ( decompressFile ( "tempFile", "tempfile2" ) );
  assert( identicalFiles("tempfile2", "tests/test2.orig"));

  assert ( compressFile ( "tests/test3.orig", "tempfile" ) );
  assert ( decompressFile ( "tempFile", "tempfile2" ) );
  assert( identicalFiles("tempfile2", "tests/test3.orig"));

  assert ( compressFile ( "tests/test4.orig", "tempfile" ) );
  assert ( decompressFile ( "tempFile", "tempfile2" ) );
  assert( identicalFiles("tempfile2", "tests/test4.orig"));*/

  /* assert(compressFile("tests/extra0.orig", "tempfile"));
   assert(decompressFile("tempfile", "tempfile2"));
   assert( identicalFiles("tempfile2", "tests/extra0.orig"));*/
  

  /*assert(compressFile("130LinesOfFileExtra7.txt", "tempfile"));
  assert(decompressFile("tempfile", "tempfile2"));
  assert( identicalFiles("tempfile2", "130LinesOfFileExtra7.txt"));

 assert(compressFile("140LinesOfFileExtra7.txt", "tempfile"));
  assert(decompressFile("tempfile", "tempfile2"));
  assert( identicalFiles("tempfile2", "140LinesOfFileExtra7.txt"));*/

  assert(compressFile("tests/extra9.orig", "tempFile"));
  assert(decompressFile("tempFile", "tempfile2"));
  assert( identicalFiles("tempfile2", "tests/extra9.orig"));

 /*assert ( compressFile ( "tests/test2.orig", "tempfile" ) );
 assert ( decompressFile ( "tempFile", "tempfile2" ) ); 
  assert( identicalFiles("tempfile2", "tests/test2.orig"));

 assert ( compressFile ( "tests/test3.orig", "tempfile" ) );
  assert ( decompressFile ( "tempFile", "tempfile2" ) ); 
  assert( identicalFiles("tempfile2", "tests/test3.orig"));

  assert ( compressFile ( "tests/test4.orig", "tempfile" ) );
  assert ( decompressFile ( "tempFile", "tempfile2" ) ); 
  assert( identicalFiles("tempfile2", "tests/test4.orig"));


  assert ( compressFile ( "tests/extra0.orig", "tempfile" ) );
  assert ( decompressFile ( "tempFile", "tempfile2" ) ); 
  assert( identicalFiles("tempfile2", "tests/extra0.orig"));

  assert ( compressFile ( "tests/extra1.orig", "tempfile" ) );
  assert ( decompressFile ( "tempFile", "tempfile2" ) ); 
  assert( identicalFiles("tempfile2", "tests/extra1.orig"));

assert ( compressFile ( "tests/extra2.orig", "tempfile" ) );
  assert ( decompressFile ( "tempFile", "tempfile2" ) ); 
  assert( identicalFiles("tempfile2", "tests/extra2.orig"));


  assert ( compressFile ( "tests/extra3.orig", "tempfile" ) );
  assert ( decompressFile ( "tempFile", "tempfile2" ) ); 
  assert( identicalFiles("tempfile2", "tests/extra3.orig"));*/


  return 0;
}
#endif /* __PROGTEST__ */
