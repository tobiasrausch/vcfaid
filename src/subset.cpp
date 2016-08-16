/*
============================================================================
VCFaid
============================================================================
Copyright (C) 2016 Tobias Rausch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Tobias Rausch (rausch@embl.de)
============================================================================
*/

#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <fstream>

#define BOOST_DISABLE_ASSERTS
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "gq.h"

using namespace vcfaid;

struct Config {
  bool hasIdFile;
  bool hasPosFile;
  boost::filesystem::path idscorefile;
  boost::filesystem::path posfile;
  boost::filesystem::path outfile;
  boost::filesystem::path vcffile;
};

template<typename TConfig, typename TScores>
inline bool
_parseScores(TConfig const& c, TScores& scores)
{
  typedef typename TScores::mapped_type TMappedType;
  typedef typename TScores::key_type TKeyType;
  std::ifstream bedFile(c.idscorefile.string().c_str(), std::ifstream::in);
  bool scoresPresent = true;
  if (bedFile.is_open()) {
    while (bedFile.good()) {
      std::string line;
      getline(bedFile, line);
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep(" \t,;");
      Tokenizer tokens(line, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter!=tokens.end()) {
	TKeyType id = *tokIter++;
	if (tokIter!=tokens.end()) {
	  TMappedType score = boost::lexical_cast<TMappedType>(*tokIter++);
	  scores[id] = score;
	} else {
	  scoresPresent = false;
	  scores[id] = 0;
	}
      }
    }
  }
  return scoresPresent;
}

template<typename TConfig, typename TGenomicPos>
inline void
_parsePositions(TConfig const& c, TGenomicPos& svpos)
{
  typedef typename TGenomicPos::value_type TChrPosPair;
  typedef typename TChrPosPair::value_type TPairSet;

  // Open VCF file
  htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);

  // Get number of sequences
  const char** seqnames = NULL;
  int32_t nseq=0;
  seqnames = bcf_hdr_seqnames(hdr, &nseq);
  if (seqnames!=NULL) free(seqnames);
  
  // Resize vectors
  svpos.resize(nseq, TChrPosPair());
  for(int32_t i = 0; i<nseq; ++i) svpos[i].resize(nseq, TPairSet());

  // Fill positions
  std::ifstream bedFile(c.posfile.string().c_str(), std::ifstream::in);
  if (bedFile.is_open()) {
    while (bedFile.good()) {
      std::string line;
      getline(bedFile, line);
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep(" \t,;");
      Tokenizer tokens(line, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter!=tokens.end()) {
	std::string chr = *tokIter++;
	if (tokIter!=tokens.end()) {
	  int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	  if (tokIter!=tokens.end()) {
	    std::string chr2 = *tokIter++;
	    if (tokIter!=tokens.end()) {
	      int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	      int32_t tid = bcf_hdr_name2id(hdr, chr.c_str());
	      int32_t mid = bcf_hdr_name2id(hdr, chr2.c_str());
	      if ((tid>=0) && (mid>=0)) svpos[tid][mid].insert(std::make_pair(start, end));
	    }
	  }
	}
      }
    }
  }

  // Close VCF
  bcf_hdr_destroy(hdr);
  bcf_close(ifile);
}

template<typename TConfig, typename TGenomicPos, typename TScores>
inline int32_t 
_processVCF(TConfig const& c, TGenomicPos const& svpos, TScores const& scores, bool hasScores) {

  // Open VCF file
  htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);

  // Open output file
  htsFile *fp = hts_open(c.outfile.string().c_str(), "wb");
  bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);
  if (hasScores) { 
    bcf_hdr_remove(hdr_out, BCF_HL_INFO, "SCORE");
    bcf_hdr_append(hdr_out, "##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"Structural Variant Score.\">");
  }
  bcf_hdr_write(fp, hdr_out);

  // Variables
  int32_t nchr2 = 0;
  char* chr2 = NULL;
  int32_t nsvend = 0;
  int32_t* svend = NULL;

  // Process records
  bcf1_t* rec = bcf_init();
  while (bcf_read(ifile, hdr, rec) == 0) {
    bcf_unpack(rec, BCF_UN_INFO);
    if (c.hasIdFile) {
      std::string svid(rec->d.id);
      if (scores.find(svid) != scores.end()) {
	if (hasScores) {
	  float score = scores.at(svid);
	  _remove_info_tag(hdr_out, rec, "SCORE");
	  bcf_update_info_float(hdr_out, rec, "SCORE", &score, 1);
	}
	
	// Write record
	bcf_write1(fp, hdr_out, rec);
      }
    } else if (c.hasPosFile) {
      bcf_get_info_string(hdr, rec, "CHR2", &chr2, &nchr2);
      std::string chr2Name(chr2);
      int32_t mid = bcf_hdr_name2id(hdr, chr2Name.c_str());
      if (mid>=0) {
	bcf_get_info_int32(hdr, rec, "END", &svend, &nsvend);
	if (svpos[rec->rid][mid].find(std::make_pair(rec->pos + 1, *svend)) != svpos[rec->rid][mid].end()) {
	  // Write record
	  bcf_write1(fp, hdr_out, rec);
	}
      }
    }
  }
  bcf_destroy(rec);
  
  // Clean-up
  if (svend != NULL) free(svend);
  if (chr2 != NULL) free(chr2);

  // Close output VCF
  bcf_hdr_destroy(hdr_out);
  hts_close(fp);
  
  // Build index
  bcf_index_build(c.outfile.string().c_str(), 14);

  // Close VCF
  bcf_hdr_destroy(hdr);
  bcf_close(ifile);
  return 0;
}

int main(int argc, char **argv) {

#ifdef PROFILE
  ProfilerStart("subset.prof");
#endif

  Config c;

  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("tsv,t", boost::program_options::value<boost::filesystem::path>(&c.idscorefile), "tab-delimited file of id & score of variants to keep")
    ("pos,p", boost::program_options::value<boost::filesystem::path>(&c.posfile), "tab-delimited file of chr, start, chr2, end of variants to keep")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("var.bcf"), "BCF output file")
    ;

  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input VCF/BCF file")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);

  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] <input.vcf.gz>" << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  } 
  
  // Check VCF file
  if (!(boost::filesystem::exists(c.vcffile) && boost::filesystem::is_regular_file(c.vcffile) && boost::filesystem::file_size(c.vcffile))) {
    std::cerr << "Input VCF/BCF file is missing " << c.vcffile.string() << std::endl;
    return 1;
  }

  // Check VCF file
  c.hasIdFile = false;
  c.hasPosFile = false;
  if (vm.count("tsv")) {
    if (!(boost::filesystem::exists(c.idscorefile) && boost::filesystem::is_regular_file(c.idscorefile) && boost::filesystem::file_size(c.idscorefile))) {
      std::cerr << "Input Identifier & Score file is missing " << c.idscorefile.string() << std::endl;
      return 1;
    }
    c.hasIdFile = true;
  } else if (vm.count("pos")) {
    if (!(boost::filesystem::exists(c.posfile) && boost::filesystem::is_regular_file(c.posfile) && boost::filesystem::file_size(c.posfile))) {
      std::cerr << "Input position file is missing " << c.posfile.string() << std::endl;
      return 1;
    }
    c.hasPosFile = true;
  } else {
    std::cerr << "Either a file listing SV identifiers or a file listing SV positions need to be specified." << std::endl;
    return 1;
  }

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Parse selected Ids and Scores or positions
  typedef std::map<std::string, double> TScores;
  typedef std::pair<int32_t, int32_t> TIntPair;
  typedef std::set<TIntPair> TPairSet;
  typedef std::vector<TPairSet> TChrPairPos;
  typedef std::vector<TChrPairPos> TGenomicPos;
  TGenomicPos svpos;
  TScores scores;
  bool hasScores = false;
  if (c.hasIdFile) hasScores = _parseScores(c, scores);
  else _parsePositions(c, svpos);

  // Filter Ids and add scores
  int r=_processVCF(c, svpos, scores, hasScores);

  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;



#ifdef PROFILE
  ProfilerStop();
#endif

  return r;
}
