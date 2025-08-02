# =========================================================================
# Git configuration
# =========================================================================
# Check where the code originates
IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  EXECUTE_PROCESS(COMMAND git ls-remote --get-url WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE GIT_ORIGIN OUTPUT_STRIP_TRAILING_WHITESPACE)
  MESSAGE(STATUS "Checking git origin: " ${GIT_ORIGIN})

  # Setup git hooks
  SET(PRECOMMIT_FILE ".githooks/pre-commit")
ENDIF()

# Perform checks only if origin points to IAG, other origins can't commit
IF("${GIT_ORIGIN}" MATCHES ".iag.uni-stuttgart.de")
  # Check if the pre-commit hooks exits
  IF (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git/hooks/pre-commit)
    # Invalid symlinks do not exist for CMake
    # > Valid inside docker images, chroot, etc.

    # Check if the user home directory exists
    IF(    DEFINED ENV{HOME})
      SET(USER_HOME_DIR $ENV{HOME})
    ELSEIF(DEFINED ENV{USERPROFILE})
      SET(USER_HOME_DIR $ENV{USERPROFILE})
    ELSE()
      SET(USER_HOME_DIR "")
    ENDIF()

    # Only attempt to create a symlink if the user has a home
    IF(USER_HOME_DIR AND NOT "${USER_HOME_DIR}" STREQUAL "/" AND IS_DIRECTORY "${USER_HOME_DIR}")
      # Create otherwise
      FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git/hooks)
      EXECUTE_PROCESS(COMMAND ln -s ${CMAKE_CURRENT_SOURCE_DIR}/${PRECOMMIT_FILE} ${CMAKE_CURRENT_SOURCE_DIR}/.git/hooks/pre-commit)
    ENDIF()
  ELSE()
    # Check if the hook is the correct symlink and warn otherwise
    EXECUTE_PROCESS(COMMAND readlink ${CMAKE_CURRENT_SOURCE_DIR}/.git/hooks/pre-commit OUTPUT_VARIABLE PRECOMMIT_LINK OUTPUT_STRIP_TRAILING_WHITESPACE)
    IF (NOT ${PRECOMMIT_LINK} MATCHES "${CMAKE_CURRENT_SOURCE_DIR}/${PRECOMMIT_FILE}")
      MESSAGE (WARNING "Custom git pre-commit hook detected. Please ensure to call ${PRECOMMIT_FILE} manually.")
    ENDIF()
  ENDIF()

  # Check if the hook actually gets loaded
  EXECUTE_PROCESS(COMMAND git config --get core.hooksPath OUTPUT_VARIABLE HOOKSPATH OUTPUT_STRIP_TRAILING_WHITESPACE WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  IF (DEFINED HOOKSPATH  AND NOT "${HOOKSPATH}" STREQUAL "" AND NOT "${HOOKSPATH}" STREQUAL ".git/hooks")
    # STRING(ASCII 27 ESCAPE)
    # MESSAGE (STATUS "${ESCAPE}[34mCustom hooks path detected. Please ensure to call ${PRECOMMIT_FILE} manually.${ESCAPE}[0m")
    MESSAGE (WARNING "Custom git hooks path detected. Please ensure to call ${PRECOMMIT_FILE} manually.")
  ENDIF()
ENDIF()
