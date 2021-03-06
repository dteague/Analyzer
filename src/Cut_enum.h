#ifndef CUTS_ENUM_H_
#define CUTS_ENUM_H_

#include <string>
#include <functional>
#include <unordered_map>

enum class CUTS {
  eGen,
  eGTau,        eGTop,        eGElec,       eGMuon,       eGZ,        eGW,       eGHiggs, eGJet,
  eRVertex,     eRMuon1,      eRMuon2,      eRElec1,      eRElec2,    eRTau1,   eRTau2,
  eRJet1,       eRJet2,       eRCenJet,     eR1stJet,     eR2ndJet,   eRBJet,   eRWjet,
    eLepPair,     eMuonPair,    eElecPair,    eMixPair,     eDiJet,     eTripVeto,
    eSusyCom,     eMET,         eNuTau,       eRTrig1,      eRTrig2, 
    Final,
    First = eGen,
  Last = Final};

/* static std::unordered_map<CUTS, std::string, EnumHash> enumNames { */
/*   {CUTS::eGen, "eGen"}, */
/*   {CUTS::eGTau, "eGTau"}, {CUTS::eGTop, "eGTop"}, {CUTS::eGElec, "eGElec"}, {CUTS::eGMuon, "eGMuon"}, {CUTS::eGZ, "eGZ"}, */
/*   {CUTS::eGW, "eGW"}, {CUTS::eGHiggs, "eGHiggs"}, {CUTS::eGJet, "eGJet"},  {CUTS::eRVertex, "eRVertex"}, */
/*   {CUTS::eRMuon1, "eRMuon1"}, {CUTS::eRMuon2, "eRMuon2"}, {CUTS::eRElec1, "eRElec1"}, {CUTS::eRElec2, "eRElec2"}, */
/*   {CUTS::eRTau1, "eRTau1"},   {CUTS::eRTau2, "eRTau2"}, {CUTS::eRJet1, "eRJet1"}, {CUTS::eRJet2, "eRJet2"}, */
/*   {CUTS::eRCenJet, "eRCenJet"}, {CUTS::eR1stJet, "eR1stJet"}, {CUTS::eR2ndJet, "eR2ndJet"}, {CUTS::eRBJet, "eRBJet"}, */
/*   {CUTS::eRWjet, "eRWjet"}, {CUTS::eDiElec, "eDiElec"}, {CUTS::eDiMuon, "eDiMuon"}, {CUTS::eDiTau, "eDiTau"}, */
/*   {CUTS::eSusyCom, "eSusyCom"}, {CUTS::eMET, "eMET"}, {CUTS::eNuTau, "eNuTau"}, {CUTS::eRTrig1, "eRTrig1"}, {CUTS::eRTrig2, "eRTrig2"} */
/*   }; */



template< typename T >
class Enum {
public:
  class Iterator {
  public:
    Iterator( int value ) :  m_value( value ) { }
    T operator*( void ) const {return (T)m_value; }
    void operator++( void ) {++m_value;}
    bool operator!=( Iterator rhs ) {return m_value != rhs.m_value;}
  private:
    int m_value;
  };
};

template< typename T >
typename Enum<T>::Iterator begin( Enum<T> ) {
  return typename Enum<T>::Iterator( (int)T::First );
}

template< typename T >
typename Enum<T>::Iterator end( Enum<T> ) {
  return typename Enum<T>::Iterator( ((int)T::Last) + 1 );
}
struct EnumHash {
  template<typename T> inline typename std::enable_if<std::is_enum<T>::value, std::size_t>::type
  operator()(const T&t) const  {
    return static_cast<std::size_t>(t);
  }
};


  #endif
