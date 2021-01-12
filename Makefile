CXX = g++
CXXFLAGS = -O3

all: structural_alignment 

structural_alignment: main.o representation_phase.o calculation_phase.o preparatory_phase.o identity_score_table.o constants.o cif_parser.o file_workflow.o geometry.o biological_structures.o block_distance.o distance_matrix.o distance_score_matrix.o scoring_matrix.o penalty_decision_matrix.o direction_matrix.o traceback.o sequence_aligner.o statistic_calculator.o
	${CXX} main.o representation_phase.o calculation_phase.o preparatory_phase.o identity_score_table.o cif_parser.o file_workflow.o geometry.o biological_structures.o block_distance.o distance_matrix.o distance_score_matrix.o scoring_matrix.o penalty_decision_matrix.o direction_matrix.o sequence_aligner.o traceback.o statistic_calculator.o -o structural_alignment

main.o: main.cpp ./src/biological_structures.h preparatory_phase.h calculation_phase.h
	${CXX} ${CXXFLAGS} main.cpp ./src/biological_structures.cpp preparatory_phase.cpp calculation_phase.cpp -c

representation_phase.o: representation_phase.h calculation_phase.h ./src/file_workflow.h ./src/biological_structures.h 
	${CXX} ${CXXFLAGS} representation_phase.cpp calculation_phase.cpp ./src/file_workflow.cpp ./src/biological_structures.cpp -c

calculation_phase.o: calculation_phase.h preparatory_phase.h ./src/block_distance.h ./src/statistic_calculator.h ./src/distance_matrix.h ./src/scoring_matrix.h
	${CXX} ${CXXFLAGS} calculation_phase.cpp preparatory_phase.cpp ./src/block_distance.cpp ./src/statistic_calculator.cpp ./src/distance_matrix.cpp ./src/scoring_matrix.cpp -c

preparatory_phase.o: preparatory_phase.h ./src/file_workflow.h ./src/cif_parser.h
	${CXX} ${CXXFLAGS} preparatory_phase.cpp ./src/constants.hpp ./src/file_workflow.cpp ./src/cif_parser.cpp ./src/biological_structures.cpp -c

identity_score_table.o: ./src/identity_score_table.cpp ./src/identity_score_table.h
	${CXX} ${CXXFLAGS} ./src/identity_score_table.cpp -c

statistic_calculator.o: ./src/statistic_calculator.cpp ./src/statistic_calculator.h
	${CXX} ${CXXFLAGS} ./src/statistic_calculator.cpp -c

sequence_aligner.o: ./src/sequence_aligner.cpp ./src/sequence_aligner.h 
	${CXX} ${CXXFLAGS} ./src/sequence_aligner.cpp -c

traceback.o: ./src/traceback.cpp ./src/traceback.h 
	${CXX} ${CXXFLAGS} ./src/traceback.cpp -c

direction_matrix.o: ./src/direction_matrix.cpp ./src/direction_matrix.h
	${CXX} ${CXXFLAGS} ./src/direction_matrix.cpp -c

penalty_decision_matrix.o: ./src/penalty_decision_matrix.cpp ./src/penalty_decision_matrix.h
	${CXX} ${CXXFLAGS} ./src/penalty_decision_matrix.cpp -c

scoring_matrix.o: ./src/scoring_matrix.cpp ./src/scoring_matrix.h
	${CXX} ${CXXFLAGS} ./src/scoring_matrix.cpp -c

distance_score_matrix.o: ./src/distance_score_matrix.cpp ./src/distance_score_matrix.h
	${CXX} ${CXXFLAGS} ./src/distance_score_matrix.cpp -c

distance_matrix.o: ./src/distance_matrix.cpp ./src/distance_matrix.h
	${CXX} ${CXXFLAGS} ./src/distance_matrix.cpp -c

block_distance.o: ./src/block_distance.cpp ./src/block_distance.h
	${CXX} ${CXXFLAGS} ./src/block_distance.cpp -c

biological_structures.o: ./src/biological_structures.cpp ./src/biological_structures.h
	${CXX} ${CXXFLAGS} ./src/biological_structures.cpp -c

geometry.o: ./src/geometry.cpp ./src/geometry.h
	${CXX} ${CXXFLAGS} ./src/geometry.cpp -c

file_workflow.o: ./src/file_workflow.cpp ./src/file_workflow.h
	${CXX} ${CXXFLAGS} ./src/file_workflow.cpp -c 

cif_parser.o: ./src/cif_parser.cpp ./src/cif_parser.h
	${CXX} ${CXXFLAGS} ./src/cif_parser.cpp -c

constants.o: ./src/constants.hpp 
	${CXX} ${CXXFLAGS} ./src/constants.hpp -c

clean:
	rm -f ./structural_alignment ./*.o





