/***********************************************************************************

    Copyright (C) 2013 by Sajia Akhter, Edwards Lab, San Diego State University

    This file is part of Qudaich.

    Qudaich is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

**********************************************************************************/
//sequence_processing.h

char codonList[4][4][4];

void init_codonList(void) {
  codonList[converter['t']][converter['t']][converter['t']] = 'F'; 
  codonList[converter['t']][converter['t']][converter['c']] = 'F';
  codonList[converter['t']][converter['t']][converter['a']] = 'L';
  codonList[converter['t']][converter['t']][converter['g']] = 'L';
  codonList[converter['c']][converter['t']][converter['t']] = 'L';
  codonList[converter['c']][converter['t']][converter['c']] = 'L';
  codonList[converter['c']][converter['t']][converter['a']] = 'L';
  codonList[converter['c']][converter['t']][converter['g']] = 'L';
  codonList[converter['a']][converter['t']][converter['t']] = 'I';
  codonList[converter['a']][converter['t']][converter['c']] = 'I';
  codonList[converter['a']][converter['t']][converter['a']] = 'I';
  codonList[converter['a']][converter['t']][converter['g']] = 'M';
  codonList[converter['g']][converter['t']][converter['t']] = 'V';
  codonList[converter['g']][converter['t']][converter['c']] = 'V';
  codonList[converter['g']][converter['t']][converter['a']] = 'V';
  codonList[converter['g']][converter['t']][converter['g']] = 'V';
  codonList[converter['t']][converter['c']][converter['t']] = 'S';
  codonList[converter['t']][converter['c']][converter['c']] = 'S';
  codonList[converter['t']][converter['c']][converter['a']] = 'S';
  codonList[converter['t']][converter['c']][converter['g']] = 'S';
  codonList[converter['c']][converter['c']][converter['t']] = 'P';
  codonList[converter['c']][converter['c']][converter['c']] = 'P';
  codonList[converter['c']][converter['c']][converter['a']] = 'P';
  codonList[converter['c']][converter['c']][converter['g']] = 'P';
  codonList[converter['a']][converter['c']][converter['t']] = 'T';
  codonList[converter['a']][converter['c']][converter['c']] = 'T';
  codonList[converter['a']][converter['c']][converter['a']] = 'T';
  codonList[converter['a']][converter['c']][converter['g']] = 'T';  
  codonList[converter['g']][converter['c']][converter['t']] = 'A';
  codonList[converter['g']][converter['c']][converter['c']] = 'A';
  codonList[converter['g']][converter['c']][converter['a']] = 'A';
  codonList[converter['g']][converter['c']][converter['g']] = 'A';
  codonList[converter['t']][converter['a']][converter['t']] = 'Y';
  codonList[converter['t']][converter['a']][converter['c']] = 'Y';
  codonList[converter['t']][converter['a']][converter['a']] = '*';
  codonList[converter['t']][converter['a']][converter['g']] = '*';
  codonList[converter['c']][converter['a']][converter['t']] = 'H';
  codonList[converter['c']][converter['a']][converter['c']] = 'H';
  codonList[converter['c']][converter['a']][converter['a']] = 'Q';
  codonList[converter['c']][converter['a']][converter['g']] = 'Q';
  codonList[converter['a']][converter['a']][converter['t']] = 'N';
  codonList[converter['a']][converter['a']][converter['c']] = 'N';
  codonList[converter['a']][converter['a']][converter['a']] = 'K';
  codonList[converter['a']][converter['a']][converter['g']] = 'K';
  codonList[converter['g']][converter['a']][converter['t']] = 'D';
  codonList[converter['g']][converter['a']][converter['c']] = 'D';
  codonList[converter['g']][converter['a']][converter['a']] = 'E';
  codonList[converter['g']][converter['a']][converter['g']] = 'E';
  codonList[converter['t']][converter['g']][converter['t']] = 'C';
  codonList[converter['t']][converter['g']][converter['c']] = 'C';
  codonList[converter['t']][converter['g']][converter['a']] = '*';
  codonList[converter['t']][converter['g']][converter['g']] = 'W';
  codonList[converter['c']][converter['g']][converter['t']] = 'R';
  codonList[converter['c']][converter['g']][converter['c']] = 'R';
  codonList[converter['c']][converter['g']][converter['a']] = 'R';
  codonList[converter['c']][converter['g']][converter['g']] = 'R';
  codonList[converter['a']][converter['g']][converter['t']] = 'S';
  codonList[converter['a']][converter['g']][converter['c']] = 'S';
  codonList[converter['a']][converter['g']][converter['a']] = 'R';
  codonList[converter['a']][converter['g']][converter['g']] = 'R';
  codonList[converter['g']][converter['g']][converter['t']] = 'G';
  codonList[converter['g']][converter['g']][converter['c']] = 'G';
  codonList[converter['g']][converter['g']][converter['a']] = 'G';
  codonList[converter['g']][converter['g']][converter['g']] = 'G';

  return;
}

