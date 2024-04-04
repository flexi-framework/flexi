SET(classes
  visuReader)

SET(headers
  ../../plugin_visu.h
  ../../pluginTypes_visu.h)

SET(private_headers
  visuReader.h)

VTK_MODULE_ADD_MODULE(visuReaderModule
  CLASSES ${classes}
  HEADERS ${headers}
  PRIVATE_HEADERS ${private_headers})

VTK_MODULE_LINK(
  visuReaderModule PUBLIC ${TOOL_NAME}lib libflexishared)

PARAVIEW_ADD_SERVER_MANAGER_XMLS(
  XMLS  visuReader.xml)
