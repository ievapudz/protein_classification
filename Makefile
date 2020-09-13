CXX = g++

all: structural_alignment

structural_alignment: main.o cif_parser.o file_workflow.o geometry.o biological_structures.o block_distance.o distance_matrix.o distance_score_matrix.o scoring_matrix.o penalty_decision_matrix.o direction_matrix.o traceback.o sequence_aligner.o statistic_calculator.o
	${CXX} main.o cif_parser.o file_workflow.o geometry.o biological_structures.o block_distance.o distance_matrix.o distance_score_matrix.o scoring_matrix.o penalty_decision_matrix.o direction_matrix.o sequence_aligner.o traceback.o statistic_calculator.o -o structural_alignment

main.o: main.cpp ./src/biological_structures.h
	${CXX} main.cpp ./src/biological_structures.cpp -c

statistic_calculator.o: ./src/statistic_calculator.cpp ./src/statistic_calculator.h
	${CXX} ./src/statistic_calculator.cpp -c

sequence_aligner.o: ./src/sequence_aligner.cpp ./src/sequence_aligner.h 
	${CXX} ./src/sequence_aligner.cpp -c

traceback.o: ./src/traceback.cpp ./src/traceback.h 
	${CXX} ./src/traceback.cpp -c

direction_matrix.o: ./src/direction_matrix.cpp ./src/direction_matrix.h
	${CXX} ./src/direction_matrix.cpp -c

penalty_decision_matrix.o: ./src/penalty_decision_matrix.cpp ./src/penalty_decision_matrix.h
	${CXX} ./src/penalty_decision_matrix.cpp -c

scoring_matrix.o: ./src/scoring_matrix.cpp ./src/scoring_matrix.h
	${CXX} ./src/scoring_matrix.cpp -c

distance_score_matrix.o: ./src/distance_score_matrix.cpp ./src/distance_score_matrix.h
	${CXX} ./src/distance_score_matrix.cpp -c

distance_matrix.o: ./src/distance_matrix.cpp ./src/distance_matrix.h
	${CXX} ./src/distance_matrix.cpp -c

block_distance.o: ./src/block_distance.cpp ./src/block_distance.h
	${CXX} ./src/block_distance.cpp -c

biological_structures.o: ./src/biological_structures.cpp ./src/biological_structures.h
	${CXX} ./src/biological_structures.cpp -c

geometry.o: ./src/geometry.cpp ./src/geometry.h
	${CXX} ./src/geometry.cpp -c

file_workflow.o: ./src/file_workflow.cpp ./src/file_workflow.h
	${CXX} ./src/file_workflow.cpp -c 

cif_parser.o: ./src/cif_parser.cpp ./src/cif_parser.h
	${CXX} ./src/cif_parser.cpp -c

clean:
	rm -f ./structural_alignment ./*.o





