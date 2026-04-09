/**
# Helper functions for output_xdmf.h
*/

#ifndef OUTPUT_XDMF_HELPERS_H
#define OUTPUT_XDMF_HELPERS_H

#include "acastillo/output_fields/output_common_helpers.h"

#ifdef HAVE_HDF5
#include "output_xdmf_helpers_xml.h"
#include "output_xdmf_helpers_populate.h"
#endif // HAVE_HDF5

#endif // OUTPUT_XDMF_HELPERS_H
