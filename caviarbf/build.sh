cd ./src
g++ -I../tclap-1.2.1/include model_search.cpp bf_io.h bf_io.cpp InputArgument.h InputArgument.cpp CombinationTable.h CombinationTable.cpp model_selection.h model_selection.cpp -o ../model_search
g++ -I../eigen-eigen-1306d75b4a21/ -I../tclap-1.2.1/include/ InputArgument.h InputArgument.cpp caviarbf.cpp -o ../caviarbf
