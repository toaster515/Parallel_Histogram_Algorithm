
include(FindPackageHandleStandardArgs)

set(charmc_compiler_names charmc)
set(charmxi_compiler_names charmxi)
set(ampi_cc_names ampicc)
set(ampi_cxx_names ampicxx)


set(possible_charm_installations ~/usr /usr/local /usr /usr/charm* /usr/local/charm* /opt/charm*)
file(GLOB possible_charm_installations ${possible_charm_installations})
#message("possible locations" ${possible_charm_installations})

#CHARM Compiler
find_program(CHARM_COMPILER
	NAMES ${charmc_compiler_names}
	HINTS ${CHARM_PATH} ${CHARM_HOME} ENV CHARM_PATH ENV CHARM_HOME
	PATHS ${possible_charm_installations}
	PATH_SUFFIXES bin
	DOC "Charm++ compiler wrapper"
)
mark_as_advanced(CHARM_COMPILER)

#Get the version
if(CHARM_COMPILER)
	find_package(MPI)

	execute_process(COMMAND ${CHARM_COMPILER} -V
				OUTPUT_VARIABLE charmc_version
				ERROR_QUIET
				OUTPUT_STRIP_TRAILING_WHITESPACE
			)
	if(charmc_version MATCHES "^Charm\\+\\+ Version [0-9.]+")
		string(REGEX REPLACE "Charm\\+\\+ Version ([0-9.]+).*" "\\1" CHARM_VERSION_STRING "${charmc_version}")
	endif()
	unset(charmc_version)
endif(CHARM_COMPILER)

#CHARMXI charm module compiler
find_program(CHARMXI_COMPILER
	NAMES ${charmxi_compiler_names}
	HINTS ${CHARM_PATH} ${CHARM_HOME} ENV CHARM_PATH ENV CHARM_HOME
	PATHS ${possible_charm_installations}
	PATH_SUFFIXES bin
	DOC "Charm++ module compiler"
)
mark_as_advanced(CHARMXI_COMPILER)

#Get all options linking, etc.
if(CHARM_COMPILER)
	execute_process(COMMAND ${CHARM_COMPILER} -print-building-blocks
				OUTPUT_VARIABLE charmc_all_variables
				ERROR_QUIET
				OUTPUT_STRIP_TRAILING_WHITESPACE
			)


	string(REPLACE "\n" ";" charmc_all_variables_list ${charmc_all_variables})
	#TODO: loop and find all such variables, not hardcoded
	##if(charmc_all_variables MATCHES "CHARM_CC_FLAGS='.*'")
	##	message(FATAL_ERROR "Found charm CC flags " ${charmc_all_variables})
	##	##string(REGEX REPLACE "Charm\\+\\+ Version ([0-9.]+).*" "\\1" CHARM_VERSION_STRING "${charmc_version}")
	##endif()
	foreach(one_charm_variable_line ${charmc_all_variables_list})
		string(REGEX REPLACE "^(.*)='(.*)'$" "\\1" ONE_CHARM_VAR_NAME ${one_charm_variable_line})
		string(REGEX REPLACE "^(.*)='(.*)'$" "\\2" ONE_CHARM_VAR_VALUE ${one_charm_variable_line})
		#message("BALAHA\n" ${ONE_CHARM_VAR_NAME} " IS EQUAL TO " ${ONE_CHARM_VAR_VALUE})
		set(${ONE_CHARM_VAR_NAME} ${ONE_CHARM_VAR_VALUE})
	endforeach()
	unset(charmc_all_variables_list)
	unset(charmc_all_variables)
	#link_directories(${target_name} ${CHARMLIB} ${CHARMLIBSO}) #Don't like that this will make everyting have that globally.
endif(CHARM_COMPILER)

if(CHARMXI_COMPILER)
	define_property(TARGET PROPERTY "CHARM_SOURCES"
		BRIEF_DOCS "Sources for charmxi"
		FULL_DOCS  "List of source files that the charm module compiler should interpret."
	)

	function(set_charm_target target_name)
		set(options SEARCH STANDALONE NOMAIN) #Tells if we want to search for .ci files in the basic sources list
		set(oneValueArgs TRACEMODE) #TODO: actually look at charmc to figure out how to properly build traces
		set(multiValueArgs CHARM_SOURCES CHARM_MODULES )
		cmake_parse_arguments(SET_CHARM_TARGET "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

		#examine all the sources and find any charm sources.
		#print_target_properties(${target_name})

		get_target_property(ALL_SOURCES_PATHS ${target_name} SOURCES)

		set(TMP_CHARM_SOURCES ) #TODO: Start with / append any charm sources provided here.
		set(TMP_NON_CHARM_SOURCES )

		foreach(one_source ${ALL_SOURCES_PATHS})
			#message("blah : " ${one_source})
			if(${one_source} MATCHES "\\.ci$")
				#message("Appending : " ${one_source})
				list(APPEND TMP_CHARM_SOURCES ${one_source})
			else()
				list(APPEND TMP_NON_CHARM_SOURCES ${one_source})
			endif()
		endforeach(one_source)
		foreach(one_charm_source ${TMP_CHARM_SOURCES})
			list(REMOVE_ITEM ALL_SOURCES_PATHS ${one_charm_source})
		endforeach(one_charm_source)

		#set_target_properties(${target_name} PROPERTIES SOURCES "${TMP_NON_CHARM_SOURCES}" SCOPE PARENT_SCOPE)
		#TODO: append to if the charm sources property already exists
		#set_target_properties(${target_name} PROPERTIES "CHARM_SOURCES" "${TMP_CHARM_SOURCES}" SCOPE PARENT_SCOPE)

		#message("all charm sources : " "${TMP_CHARM_SOURCES}")
		#message("all non-charm sources : " "${TMP_NON_CHARM_SOURCES}")

		#message("THE CHARMXI COMPILER IS " ${CHARMXI_COMPILER})

		foreach(one_charm_source ${TMP_CHARM_SOURCES})
			get_filename_component(SINGLE_CHARM_DEFAULT_OUTPUT ${one_charm_source} NAME)
			string(REGEX REPLACE "\\.ci$" "" SINGLE_CHARM_DEFAULT_OUTPUT ${SINGLE_CHARM_DEFAULT_OUTPUT})

			#TODO: We should create a directory that these generated files go into.
			#If only certain modules were asked for, we should generate those into a non-default directory.

			list(APPEND TMP_NON_CHARM_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/${SINGLE_CHARM_DEFAULT_OUTPUT}.decl.h")
			list(APPEND TMP_NON_CHARM_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/${SINGLE_CHARM_DEFAULT_OUTPUT}.def.h")
			include_directories(${target_name} ${CMAKE_CURRENT_BINARY_DIR})

			#If we use an OUTPUT type custom_command, and alter the target's sources list, we might avoid that.
			#message("one_charm_source : " ${CMAKE_CURRENT_SOURCE_DIR}/${one_charm_source} )
			set(SET_CHARM_TARGET_SINGLE_CHARM_SOURCE_FULL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/${one_charm_source})
			add_custom_command(
				PRE_BUILD
				OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${SINGLE_CHARM_DEFAULT_OUTPUT}.decl.h ${CMAKE_CURRENT_BINARY_DIR}/${SINGLE_CHARM_DEFAULT_OUTPUT}.def.h
				COMMAND ${CHARMXI_COMPILER} ${SET_CHARM_TARGET_SINGLE_CHARM_SOURCE_FULL_PATH}
				WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
				DEPENDS ${SET_CHARM_TARGET_SINGLE_CHARM_SOURCE_FULL_PATH}
				VERBATIM
			)

			#TODO: needs to be per module
			set(mod_init_src "${CMAKE_CURRENT_BINARY_DIR}/${SINGLE_CHARM_DEFAULT_OUTPUT}_modinit.C")
			list(APPEND TMP_NON_CHARM_SOURCES ${mod_init_src})
			add_custom_command(
				PRE_BUILD
				OUTPUT ${mod_init_src}
				COMMAND echo "void _registerExternalModules(char **argv) { (void)argv; } void _createTraces(char **argv) {(void)argv;}" >> ${mod_init_src}
				WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
				VERBATIM
			)


			#todo: each module needs a modInit
		endforeach()

		set_target_properties(${target_name} PROPERTIES SOURCES "${TMP_NON_CHARM_SOURCES}" SCOPE PARENT_SCOPE)
		#TODO: append to if the charm sources property already exists
		set_target_properties(${target_name} PROPERTIES "CHARM_SOURCES" "${TMP_CHARM_SOURCES}" SCOPE PARENT_SCOPE)

		#compile and linking flags
		#TODO: Charm can be built without MPI, can't it?
		#TODO: Detect the language C/CXX etc.
		#TODO: Get the last compiler/linker flags dynamically from interrogating charmc, not hardcoded as they are here "-m64 etc."
		include_directories(${target_name} ${MPI_CXX_INCLUDE_PATH})
		target_link_libraries(${target_name} ${MPI_CXX_LIBRARIES})
		set_target_properties(${target_name} PROPERTIES COMPILE_FLAGS "${CHARM_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS} -m64 -fPIC " SCOPE PARENT_SCOPE)
		set_target_properties(${target_name} PROPERTIES LINK_FLAGS "${CHARM_LDXX_FLAGS} ${MPI_CXX_LINK_FLAGS} -m64 -fPIC -rdynamic " SCOPE PARENT_SCOPE)

		include_directories(${target_name} ${CHARMINC})

	endfunction()
endif(CHARMXI_COMPILER)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Charm
	FOUND_VAR "CHARM_FOUND"
	REQUIRED_VARS CHARM_COMPILER CHARMXI_COMPILER
	VERSION_VAR CHARM_VERSION_STRING)

#Also find AMPI

#AMPI C Compiler
find_program(AMPI_C_COMPILER
	NAMES ${ampi_cc_names}
	HINTS ${CHARM_PATH} ${CHARM_HOME} ENV CHARM_PATH ENV CHARM_HOME
	PATHS ${possible_charm_installations}
	PATH_SUFFIXES bin
	DOC "AMPI C compiler wrapper"
)
mark_as_advanced(AMPI_C_COMPILER)

#AMPI CXX Compiler
find_program(AMPI_CXX_COMPILER
	NAMES ${ampi_cxx_names}
	HINTS ${CHARM_PATH} ${CHARM_HOME} ENV CHARM_PATH ENV CHARM_HOME
	PATHS ${possible_charm_installations}
	PATH_SUFFIXES bin
	DOC "AMPI CXX compiler wrapper"
)
mark_as_advanced(AMPI_CXX_COMPILER)

find_package_handle_standard_args(AMPI
	FOUND_VAR "AMPI_FOUND"
	REQUIRED_VARS AMPI_C_COMPILER AMPI_CXX_COMPILER
	VERSION_VAR CHARM_VERSION_STRING)
