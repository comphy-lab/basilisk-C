export BASILISK := $(patsubst %/src/build-env/Makefile,%,$(realpath $(lastword $(MAKEFILE_LIST))))

check: minimal gview gpus

# $(HOME)/.cache/podman/apt needs to exist, mkdir it if not
APT = $(HOME)/.cache/podman/apt:/var/cache/apt

## A minimal Debian + gcc an make

basilisk-minimal:
	podman build -v "$(APT)" -t basilisk-minimal -f minimal.Containerfile

MINIMAL = cd /root && \
	darcs get /basilisk && \
	cd basilisk/src && \
	ln -s config.gcc config && \
	(make -k || true) && \
	make && \
	./qcc --version

minimal: basilisk-minimal
	podman run --rm \
		-v $(BASILISK):/basilisk:ro \
		basilisk-minimal \
		bash -c '$(MINIMAL)'

## Graphical user interfaces

basilisk-gview:
	podman build -v "$(APT)" -t basilisk-gview -f gview.Containerfile

GVIEW = $(MINIMAL) && \
	cd gview && \
	make && \
	make gview-nogui gview

gview: basilisk-gview
	podman run --rm \
		-v $(BASILISK):/basilisk:ro \
		basilisk-gview \
		bash -c '$(GVIEW)'

## GPU grids

basilisk-gpus:
	podman build -v "$(APT)" -t basilisk-gpus -f gpus.Containerfile

GPUS = $(MINIMAL) && \
	cd grid/gpu && make libgpu.a && \
	cd ../cuda && make libbuda.a && \
	cd ../hip && make libhipcuda.a libhipamd.a && \
	cd ../opencl && make libocl.a

gpus: basilisk-gpus
	podman run --rm \
		-v $(BASILISK):/basilisk:ro \
		basilisk-gpus \
		bash -c '$(GPUS)'

# To run interactively do:
#  podman run --rm -it -v ~/basilisk:/basilisk basilisk-minimal bash
