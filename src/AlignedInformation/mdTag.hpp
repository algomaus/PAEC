
#pragma once

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>
#include <vector>
#include <string>

typedef struct MDEntry {
  seqan::Dna substitutedNucleotide;
  unsigned position;
} MDEntry;

typedef std::vector<MDEntry> MDTag;

inline std::string mdTagToString(MDTag &mdTag) {
	std::string res = "";
	for (size_t i = 0; i < mdTag.size(); ++i) {
		res += mdTag[i].substitutedNucleotide;
		res += " " + std::to_string(mdTag[i].position);
		res += ",";
	}
	return res;
}

inline void parseMDTag(seqan::CharString mdString, MDTag &mdTag) {
  bool delmode = false;
  
  //std::cout << "mdString: " << mdString << std::endl;
  
  unsigned position = 0;
  std::string tempString = "";
  for (unsigned i = 0; i < length(mdString); ++i) {
    if (mdString[i] == 'A' || mdString[i] == 'C' || mdString[i] == 'G' || mdString[i] == 'T') {
      if (delmode) {
      } else {
	if (!tempString.empty()) {
	  position += std::stoi(tempString);
	}
	MDEntry entry;
	entry.substitutedNucleotide = mdString[i];
	entry.position = position;
	mdTag.push_back(entry);
	tempString.clear();
	position++;
      }
    } else if (mdString[i] == '^') { // deletion
      if (!tempString.empty()) {
	position += std::stoi(tempString);
      }
      tempString.clear();
      delmode = true;
    } else { // number
      delmode = false;
      tempString += mdString[i];
    }
  }
}

inline void extractMDTag(seqan::BamAlignmentRecord &record, MDTag &mdTag) {
  seqan::CharString mdString;
  seqan::BamTagsDict tagsDict(record.tags);
  unsigned tagIdx = 0;
  if (!findTagKey(tagIdx, tagsDict, "MD")) {
    std::cerr << "ERROR: Unknown key!\n";
  }
  if (!extractTagValue(mdString, tagsDict, tagIdx)) {
    std::cerr << "ERROR: There was an error extracting MD from tags!\n";
  }
  
  parseMDTag(mdString, mdTag);
}
