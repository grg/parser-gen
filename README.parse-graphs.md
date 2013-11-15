Parse graph format
==================

The parse graph file describes the parse graph. It does so by listing each of
the header types and how they connect with one another. The format of the
header file is as follows:

    header1 {
        fields {
            field1 : len : extract?,
            ...
            fieldN : len : extract?,
        }
    
        pseudo-fields {
            pseudoField1 : len,
            ...
            pseudoFieldK : len,
        }
    
        next_header = NEXT_HEADER_TYPE
            OR
        next_header = map(lookupField1, ..., lookupFieldM) {
            mapValue1 : NEXT_HEADER_TYPE1,
            mapValue2 : NEXT_HEADER_TYPE2,
            ...
            mapValueK : NEXT_HEADER_TYPEK,
        }
        next_header_def = ???
    
        max_var = ???
        max = ???
        length = ???
        max_length = ???
    }
    
    header2 {
        # Header 2 description
    }
    
    ...
    
    headerN {
        # Header N description
    }

The names of the headers should appear in place of the `headerX` values.

`fields` describes each of the fields within the header. Specify the name and
length for each field, and, optionally, add the word `extract` if the field
should be extracted into the packet header vector. Fields are listed in the
order they appear in the header. Variable length headers end with a field with
a length specified as an asterisk.

`pseudo-fields` is optional and describes any fields from the _next_ header
that are used when making next-header/length decisions. Pseudo-fields are used
in decision only and are not extracted or consumed. An example use case for
pseudo-fields is the MPLS header; in pseudo-wire configurations, the four bits
after the MPLS header identify whether the next header type is IPv4 (value:
4), IPv6 (value: 6), or Ethernet-over-MPLS (value: 0).

`next_header` identifies the header type(s) the follow this header. A given
header may either be followed _always_ by one particular header or may be
followed by several possible header types determined by one or more fields
wihtin the header. In the former case, the `next_header` line specifies the
header type that always follows. In the latter case, the `next_header` line
contains the word `map` followed by a list of fields that are used to make the
decision. The fields are concatenated together and this concatenated value is
used to identify the next header type. Following the `next_header` line are a
set of values and the correspond next header types.

`next_header_def` is optional and specifies the default value to use for the
`next_header` input value to use in the test bench when not following any of
the transitions.

`max_var` and `max` are optional. They are used when a header may occur
multiple times up to a maximum count. The `max` value specifies the maximum
count; the `max_var` is a name given to the counter. Multiple headers may
share a counter to limit the combined total of those headers.

`length` is an optional field for variable-length headers. It specifies how
the length of the header is determined using the fields within the packet.

`max_length` is an optional field that specifies the _maximum_ length of a
variable-length header.

Comments are preceded by a hash (#).


Example header definitions
==========================

Ethernet
--------

The Ethernet header description is:

    ethernet {
        fields {
            dstAddr : 48 : extract,
            srcAddr : 48 : extract,
            etherType : 16 : extract,
        }
        next_header = map(etherType) {
            0x8100, 0x9100 : ieee802-1q,
            0x8847, 0x8848 : mpls,
            0x0800 : ipv4,
            0x86dd : ipv6,
            0x0806, 0x8035 : arp_rarp,
        }
    }

The Ethernet header consists of three fields:
  * dstAddr: a 48-bit field to be extracted
  * srcAddr: a 48-bit field to be extracted
  * etherType: a 16-bit field to be extracted

The next header is determined by the etherType field value. For example, a
value of 0x8100 indicates the next header type is IEEE 802.1Q (a VLAN tag),
while 0x0800 indicates the next header type is IPv4.


IPv4
----

    ipv4 {
        fields {
            version : 4,
            ihl : 4,
            diffserv : 8 : extract,
            totalLen : 16,
            identification : 16,
            flags : 3 : extract,
            fragOffset : 13,
            ttl : 8 : extract,
            protocol : 8 : extract,
            hdrChecksum : 16,
            srcAddr : 32 : extract,
            dstAddr : 32 : extract,
            options : *,
        }
        next_header = map(fragOffset, protocol) {
            1 : icmp,
            6 : tcp,
            17 : udp,
            47 : gre,
            50 : ipsec_esp,
            51 : ipsec_ah,
            132 : sctp,
        }
        length = ihl * 4 * 8
        max_length = 256
    }

The IPv4 header is variable-length as indicated by the options field with a
length of asterisk. Only a subset of the fields are extracted and placed in
the packet header vector.

The next header type is determined by a combination of the fragment offset and
protocol fields. For example, a value of 1 identifies ICMP as the next header
type. (The fragment offset is 0 for all the transitions.)

The header definition declares the length as calculated as the `ihl` field
multiplied by 4 * 8. (The IHL field indicates the number of 32-bit/4-byte
words in the header.)
