CC = gcc
CFLAGS = -Wall -O2 -fopenmp -lm
TARGET = simulador
SRC = main.c init.c utils.c

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET) snapshots/*.bin

run: $(TARGET)
	./$(TARGET)

.PHONY: clean run