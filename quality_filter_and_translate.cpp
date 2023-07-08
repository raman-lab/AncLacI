//File Name: quality_filter_and_translate.cpp
//Author: Anthony Meger
//Email Address: ameger@wisc.edu
//Description: program to analyze NGS pair-end reads; evaluates read quality, identifies barcode and performs translations
//Last Changed: 8/2/2021

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
using namespace std;

string dna2aa (string dna)
{
  for(int i=0; i < 3; i++)
  {
    if((dna.substr(i,1) != "A") & (dna.substr(i,1) != "C") & (dna.substr(i,1) != "T") & (dna.substr(i,1) != "G"))
    {
      return("?");
    } 
  } 
  if(dna == "AAA")
  {
    return("K");
  }
  if(dna == "AAC")
  {
    return("N");
  }
  if(dna == "AAG")
  {
    return("K");
  }
  if(dna == "AAT")
  {
    return("N");
  }
  if(dna == "ACA")
  {
    return("T");
  }
  if(dna == "ACC")
  {
    return("T");
  }
  if(dna == "ACG")
  {
    return("T");
  }
  if(dna == "ACT")
  {
    return("T");
  }
  if(dna == "AGA")
  {
    return("R");
  }
  if(dna == "AGC")
  {
    return("S");
  }
  if(dna == "AGG")
  {
    return("R");
  }
  if(dna == "AGT")
  {
    return("S");
  }
  if(dna == "ATA")
  {
    return("I");
  }
  if(dna == "ATC")
  {
    return("I");
  }
  if(dna == "ATG")
  {
    return("M");
  }
  if(dna == "ATT")
  {
    return("I");
  }
  if(dna == "CAA")
  {
    return("Q");
  }
  if(dna == "CAC")
  {
    return("H");
  }
  if(dna == "CAG")
  {
    return("Q");
  }
  if(dna == "CAT")
  {
    return("H");
  }
  if(dna == "CCA")
  {
    return("P");
  }
  if(dna == "CCC")
  {
    return("P");
  }
  if(dna == "CCG")
  {
    return("P");
  }
  if(dna == "CCT")
  {
    return("P");
  }
  if(dna == "CGA")
  {
    return("R");
  }
  if(dna == "CGC")
  {
    return("R");
  }
  if(dna == "CGG")
  {
    return("R");
  }
  if(dna == "CGT")
  {
    return("R");
  }
  if(dna == "CTA")
  {
    return("L");
  }
  if(dna == "CTC")
  {
    return("L");
  }
  if(dna == "CTG")
  {
    return("L");
  }
  if(dna == "CTT")
  {
    return("L");
  }
  if(dna == "GAA")
  {
    return("E");
  }
  if(dna == "GAC")
  {
    return("D");
  }
  if(dna == "GAG")
  {
    return("E");
  }
  if(dna == "GAT")
  {
    return("D");
  }
  if(dna == "GCA")
  {
    return("A");
  }
  if(dna == "GCC")
  {
    return("A");
  }
  if(dna == "GCG")
  {
    return("A");
  }
  if(dna == "GCT")
  {
    return("A");
  }
  if(dna == "GGA")
  {
    return("G");
  }
  if(dna == "GGC")
  {
    return("G");
  }
  if(dna == "GGG")
  {
    return("G");
  }
  if(dna == "GGT")
  {
    return("G");
  }
  if(dna == "GTA")
  {
    return("V");
  }
  if(dna == "GTC")
  {
    return("V");
  }
  if(dna == "GTG")
  {
    return("V");
  }
  if(dna == "GTT")
  {
    return("V");
  }
  if(dna == "TAA")
  {
    return("Z");
  }
  if(dna == "TAC")
  {
    return("Y");
  }
  if(dna == "TAG")
  {
    return("Z");
  }
  if(dna == "TAT")
  {
    return("Y");
  }
  if(dna == "TCA")
  {
    return("S");
  }
  if(dna == "TCC")
  {
    return("S");
  }
  if(dna == "TCG")
  {
    return("S");
  }
  if(dna == "TCT")
  {
    return("S");
  }
  if(dna == "TGA")
  {
    return("Z");
  }
  if(dna == "TGC")
  {
    return("C");
  }
  if(dna == "TGG")
  {
    return("W");
  }
  if(dna == "TGT")
  {
    return("C");
  }
  if(dna == "TTA")
  {
    return("L");
  }
  if(dna == "TTC")
  {
    return("F");
  }
  if(dna == "TTG")
  {
    return("L");
  }
  if(dna == "TTT")
  {
    return("F");
  }
}

int getReadingFrame (string seq, const string read_octet, int start, int end) 
{
  int start_position = 0;
  string octet;

  for(int nt = start; nt < end; nt++)
  {
    octet = seq.substr(nt,8);
    if(octet == read_octet)
    {
      start_position = nt + 8;
      return(start_position);
    } 
  }
  return(0);
}

string getTranslation (string seq, int start_pos)
{
  int aa_len = floor((seq.length() - start_pos) / 3);
  string aa_seq = "", codon;

  for(int pos=0; pos < aa_len; pos++)
  {
    codon = seq.substr((start_pos + (3 * pos)),3);
    aa_seq = aa_seq + dna2aa(codon);
    if(dna2aa(codon) == "Z")
    {
      return(aa_seq);
    }
  }
  return(aa_seq);
}

int getFileLength(string filename)
{
  int lineCount=0;
  string line;
  ifstream inFile(filename.c_str());
  while(getline(inFile,line))
  {
    lineCount++;
  }
  return(lineCount);
}

int main()
{
  cout << "This script should be used after merging paired-end reads." << '\n';
  cout << "Please ensure that your input file is formatted like this:" << '\n';
  cout << "  DNAsequence1 QSCOREsequence1" << '\n';
  cout << "  DNAsequence2 QSCOREsequence2" << '\n';
  cout << "  ...          ..." << '\n';
  
  string read_octet="AAGGTACC";
  int rf_scan_start=0, rf_scan_end;
  string inMPERs="combined.fastq";
  ifstream inRawReads(inMPERs.c_str());
  cout << "Processing file now..." << '\n';

  ofstream outGoodReads("good_reads.csv"), outPoorReads("poor_reads.csv");
  double p, q_adj, q;
  string line, str, seq;
  int start=0;

  while(getline(inRawReads,line))
  {
    for(int i=0; i < line.length(); i++)
    {
      str = line.substr(i,1);
      if(str == " ")
      {
        str = line.substr(i+1,(line.length()-i-1));
        seq = line.substr(0,line.length()-i-1);
        break;  
      }
    }
    p = 0.0;

    for(int j=0; j < str.length(); j++)
    {
      q = double(str[j]) - 33;
      q_adj = (-1 * q) / 10;
      p = p + pow(10,q_adj);
    }

    if(p < 10.0)
    {
      outGoodReads << seq << '\n';
    } else
    {
      outPoorReads << seq << "," << p << '\n';
    }
  }
  inRawReads.close();
  outGoodReads.close();
  outPoorReads.close();

  int numGoodReads=getFileLength("good_reads.csv");
  ifstream inGoodReads("good_reads.csv");
  ofstream translations("translations.csv");
  ofstream failed("failed_to_translate.csv");
  string barcode="none", dna_sequence, aa_sequence;
  int tsl_start;

  while(getline(inGoodReads, dna_sequence))
  {
    rf_scan_end = dna_sequence.length() - 8;    
    tsl_start = getReadingFrame(dna_sequence, read_octet, rf_scan_start, rf_scan_end);
    if(tsl_start > 0)
    {
      aa_sequence = getTranslation(dna_sequence, tsl_start);
      translations << aa_sequence << "," << barcode << '\n';
    } else
    {
      failed << dna_sequence << "," << barcode << '\n';
    }
  }
      
  inGoodReads.close();
  failed.close();
  translations.close();
  
  cout << "DNA sequences passing the quality filter were written to 'good_reads.csv'" << '\n';
  cout << "DNA sequences failing the quality filter were written to 'poor_read.csv'" << '\n';
  cout << "DNA Sequences lacking a defined reading frame were written to 'failed_to_translate.csv'" << '\n';
  cout << "Translated amino acid sequences were written to 'translations.csv'" << '\n';
  cout << "Note, 'z' denotes a termination codon and '?' denotes an ambiguous base call " << '\n';

  return 0;
}

