
#pragma once

CorrectedRead precorrectRead_KmerBased(const FASTQRead &fastqRead, ErrorProfileUnit &errorProfile, KmerClassificationUnit &kmerClassifier);
CorrectedRead precorrectRead_Naive(const FASTQRead &fastqRead, ErrorProfileUnit &errorProfile);
CorrectedRead postcorrectRead_Multidel(const FASTQRead &fastqRead, ErrorProfileUnit &errorProfile, KmerClassificationUnit &kmerClassifier);
