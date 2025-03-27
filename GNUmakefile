CC := gcc
CFLAGS := -Wall -Wextra -Wno-builtin-declaration-mismatch
LDFLAGS := -lm

SRC := src/materials.c src/scene.c src/compute_paths.c test/test.c viz/vizrays.c
SCENE_FROMSIONNA_SRC := src/scene.c src/materials.c src/scene_fromSionna.c

TARGET := test.elf
SCENE_FROMSIONNA_TARGET := scene_fromSionna.elf

.PHONY: all dbg clean scene_fromSionna scene_fromSionna_dbg

all: CFLAGS += -O3
all: $(TARGET)

dbg: CFLAGS += -g -O0
dbg: $(TARGET)

scene_fromSionna: CFLAGS += -O3
scene_fromSionna: $(SCENE_FROMSIONNA_TARGET)

scene_fromSionna_dbg: CFLAGS += -g -O0
scene_fromSionna_dbg: $(SCENE_FROMSIONNA_TARGET)

$(TARGET): LDFLAGS += -lGL -lGLU -lglut
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

$(SCENE_FROMSIONNA_TARGET): $(SCENE_FROMSIONNA_SRC)
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm -f $(TARGET) $(SCENE_FROMSIONNA_TARGET) rt rt.*.so *.bin *.elf *.o
	rm -rf __pycache__
