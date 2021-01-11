#ifndef _REPRESENTATION_PHASE_
#define _REPRESENTATION_PHASE_
#include "calculation_phase.h"
#include "./src/biological_structures.h"
#include "./src/file_workflow.h"

class RepresentationPhase{
        CalculationPhase* calculations_;
        Protein& p_protein_;
        Protein& q_protein_;
    public:
        RepresentationPhase(CalculationPhase* calculations);
        void representAlignment();
        void representNumeralAlignment();
        std::pair<double, double> representIdentityScore();
};


#endif
