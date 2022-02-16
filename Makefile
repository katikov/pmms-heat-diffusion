#MODIFY the following variables (STUDENTIDS, GROUP) with your infos
#STUDENTIDS should be a list of the student ids of all group members,
#			e.g, STUDENTIDS=1111111 2222222 3333333 4444444
STUDENTIDS=14113139 13676105 14100932 14122073
GROUP=11

ASSIGNMENTS=$(filter assignment_%,$(wildcard *))
EMPTY :=
SPACE := $(EMPTY) $(EMPTY)

.PHONY: $(ASSIGNMENTS)
all: $(ASSIGNMENTS)

$(ASSIGNMENTS):
	$(MAKE) -C $@

CFLAGS=-I./include
.PHONY: clean submission demo

demo:
	$(MAKE) -C demo

$(basename $(filter %.c,$(wildcard demo/*)) ): $(filter-out obj/main.o, $(OBJ))
	gcc $(CFLAGS) $@.c $(filter-out obj/main.o, $(OBJ)) -o $@ $(LIBS)

submission_%:
	$(MAKE) -C assignment_$(patsubst submission_%,%,$@) clean
	tar czfh heat_assignment_$(patsubst submission_%,%,$@)_group_$(GROUP)_$(subst $(SPACE),_,$(shell echo $(STUDENTIDS) | tr -s " ")).tar.gz assignment_$(patsubst submission_%,%,$@)

clean:
	-rm -f *.tar.gz
	$(MAKE) -C demo clean
	$(foreach assignment,$(ASSIGNMENTS),$(MAKE) -C $(assignment) clean;) 
