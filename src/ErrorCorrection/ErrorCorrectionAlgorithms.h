
#pragma once

CorrectedRead correctRead_KmerBased(const FASTQRead &fastqRead, ErrorProfileUnit &errorProfile, KmerClassificationUnit &kmerClassifier, bool correctIndels = true);
CorrectedRead correctRead_Naive(const FASTQRead &fastqRead, ErrorProfileUnit &errorProfile, KmerClassificationUnit &kmerClassifier, bool correctIndels = true);
//CorrectedRead postcorrectRead_Multidel(const FASTQRead &fastqRead, ErrorProfileUnit &errorProfile, KmerClassificationUnit &kmerClassifier);
