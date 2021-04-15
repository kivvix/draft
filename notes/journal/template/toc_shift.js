#!/usr/bin/env node

var fs = require("fs");
const jsdom = require("jsdom");
const { JSDOM } = jsdom;

const argv = require('yargs')
  .scriptName("toc_shift")
  .usage("$0 -i inputfile -o outputfile")
  .option('input',{
    alias: 'i',
    describe: 'input file where need to move `#TOC` outside of `<main>` and add `.markdown-body` class to `<main>`'
  })
  .option('output',{
    alias: 'o',
    describe: 'output file with all modifications'
  })
  .demandOption(['input','output'],"Please provide both input and output arguments to work this tool")
  .help()
  .alias('help','h')
  .argv

JSDOM.fromFile( argv['input'] ,{}).then(dom => {
  // move #TOC
  var toc  = dom.window.document.querySelector("#TOC");
  var main = dom.window.document.querySelector("main");

  main.parentNode.insertBefore(toc,main);
  main.classList.add("markdown-body");

  // write file
  fs.writeFile(argv['output'], dom.serialize(), (err) => {
    if (err) console.log(err);
    console.log("Successfully KaTeX rended!");
  });

});
