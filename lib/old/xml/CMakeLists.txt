add_definitions(-DTR15581)

add_library(xml
dom/flib_dom.f90
dom/m_dom_attribute.f90
dom/m_dom_debug.f90
dom/m_dom_document.f90
dom/m_dom_element.f90
dom/m_dom_error.f90
dom/m_dom_namednodemap.f90
dom/m_dom_node.f90
dom/m_dom_nodelist.f90
dom/m_dom_parse.f90
dom/m_dom_types.f90
dom/m_dom_utils.f90

sax/flib_sax.f90
sax/m_buffer.f90
sax/m_charset.f90
sax/m_converters.f90
sax/m_debug.f90
sax/m_dictionary.f90
sax/m_elstack.f90
sax/m_entities.f90
sax/m_fsm.f90
sax/m_io.f90
sax/m_reader.f90
sax/m_xml_error.f90
sax/m_xml_parser.f90

strings/m_strings.F90

wxml/flib_wxml.f90
wxml/m_wxml_buffer.f90
wxml/m_wxml_core.f90
wxml/m_wxml_dictionary.f90
wxml/m_wxml_elstack.f90
wxml/m_wxml_text.f90
)
