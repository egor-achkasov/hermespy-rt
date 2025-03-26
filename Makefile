CC := gcc
CFLAGS := -Wall -Wextra -Wno-builtin-declaration-mismatch
LDFLAGS := -lm

SRC := src/materials.c src/scene.c src/compute_paths.c test/test.c
VIZRAYS_SRC := src/scene.c src/vizrays.c

TARGET := test.elf
VIZRAYS_TARGET := vizrays.elf

.PHONY: all dbg clean vizrays vizraysdbg

all: CFLAGS += -O3
all: $(TARGET)

dbg: CFLAGS += -g -O0
dbg: $(TARGET)

vizrays: CFLAGS += -O3
vizrays: LDFLAGS += -lGL -lGLU -lglut
vizrays: $(VIZRAYS_TARGET)

vizraysdbg: CFLAGS += -g -O0
vizraysdbg: LDFLAGS += -lGL -lGLU -lglut
vizraysdbg: $(VIZRAYS_TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

$(VIZRAYS_TARGET): $(VIZRAYS_SRC)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

clean:
	rm -f $(TARGET) $(VIZRAYS_TARGET) rt rt.*.so *.bin a.out vizrays
	rm -rf __pycache__

