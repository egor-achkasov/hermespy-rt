CC := gcc
CFLAGS := -Wall -Wextra -Wno-builtin-declaration-mismatch
LDFLAGS := -lm

SRC := src/materials.c src/scene.c src/compute_paths.c test/test.c
VIZRAYS_SRC := src/scene.c src/vizrays.c
SCENE_FROMSIONNA_SRC := src/scene.c src/materials.c src/scene_fromSionna.c

TARGET := test.elf
VIZRAYS_TARGET := vizrays.elf
SCENE_FROMSIONNA_TARGET := scene_fromSionna.elf

.PHONY: all dbg clean vizrays vizrays_dbg scene_fromSionna scene_fromSionna_dbg

all: CFLAGS += -O3
all: $(TARGET)

dbg: CFLAGS += -g -O0
dbg: $(TARGET)

vizrays: CFLAGS += -O3
vizrays: LDFLAGS += -lGL -lGLU -lglut
vizrays: $(VIZRAYS_TARGET)

vizrays_dbg: CFLAGS += -g -O0
vizrays_dbg: LDFLAGS += -lGL -lGLU -lglut
vizrays_dbg: $(VIZRAYS_TARGET)

scene_fromSionna: CFLAGS += -O3
scene_fromSionna: $(SCENE_FROMSIONNA_TARGET)

scene_fromSionna_dbg: CFLAGS += -g -O0
scene_fromSionna_dbg: $(SCENE_FROMSIONNA_TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

$(VIZRAYS_TARGET): $(VIZRAYS_SRC)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

$(SCENE_FROMSIONNA_TARGET): $(SCENE_FROMSIONNA_SRC)
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm -f $(TARGET) $(VIZRAYS_TARGET) rt rt.*.so *.bin *.elf
	rm -rf __pycache__

