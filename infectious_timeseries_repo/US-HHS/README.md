{\rtf1\ansi\ansicpg1252\cocoartf2761
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;\f1\fswiss\fcharset0 Helvetica-Bold;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;\red24\green26\blue30;\red38\green38\blue38;
\red255\green255\blue255;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;\cssrgb\c12157\c13725\c15686;\cssrgb\c20000\c20000\c20000;
\cssrgb\c100000\c100000\c100000;}
\margl1440\margr1440\vieww15900\viewh9400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs32 \cf0 US HHS Data\

\fs24 \

\f1\b \cf2 Download Link
\f0\b0 \
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
{\field{\*\fldinst{HYPERLINK "https://github.com/reichlab/flusion/blob/main/data-raw/influenza-hhs/download-hhs.R"}}{\fldrslt \cf2 https://github.com/reichlab/flusion/blob/main/data-raw/influenza-hhs/download-hhs.R}}
\f1\b \
\
Citation\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\b0 \cf2 @misc\{ray2024flusionintegratingmultipledata,\
      title=\{Flusion: Integrating multiple data sources for accurate influenza predictions\}, \
      author=\{Evan L. Ray and Yijin Wang and Russell D. Wolfinger and Nicholas G. Reich\},\
      year=\{2024\},\
      eprint=\{2407.19054\},\
      archivePrefix=\{arXiv\},\
      primaryClass=\{stat.ML\},\
      url=\{https://arxiv.org/abs/2407.19054\}, \
\}
\f1\b \cf0 \
\
Licensing 
\f0\b0 \
\pard\pardeftab720\partightenfactor0
\cf3 \expnd0\expndtw0\kerning0
MIT License\
\
Copyright (c) 2023 The Reich Lab at UMass-Amherst\
\
Permission is hereby granted, free of charge, to any person obtaining a copy\
of this software and associated documentation files (the "Software"), to deal\
in the Software without restriction, including without limitation the rights\
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\
copies of the Software, and to permit persons to whom the Software is\
furnished to do so, subject to the following conditions:\
\
The above copyright notice and this permission notice shall be included in all\
copies or substantial portions of the Software.\
\
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\
SOFTWARE.\
\

\f1\b Data Cleaning\

\f0\b0 \
These data were obtained via the Reich Lab GitHub. Daily reported hospitalization dates were modified by 1 following precedent in the code provided by the Reich Lab. Daily reports were aggregated to weekly reports also following code provided by the Reich Lab. If a week had no reported influenza hospitalizations for a given location, we report 0 hospitalizations, ensuring that the resulting time series has regular weekly reporting without gaps. We require there to be at least 10 observations per time series for inclusion. \
\
\pard\pardeftab720\partightenfactor0

\f1\b \cf4 \cb5 Date of Raw Data Access
\f0\b0 : \
2/23/2026\cf3 \cb1 \
}