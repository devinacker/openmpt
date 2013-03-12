/*
 * PatternClipboard.h
 * ------------------
 * Purpose: Implementation of the pattern clipboard mechanism
 * Notes  : (currently none)
 * Authors: OpenMPT Devs
 * The OpenMPT source code is released under the BSD license. Read LICENSE for more details.
 */


#pragma once

#include <deque>
#include "Snd_defs.h"
#include "PatternCursor.h"

struct ModCommandPos;
class CSoundFile;

//===========================
class PatternClipboardElement
//===========================
{
public:
	
	CString content;
	CString description;

public:

	PatternClipboardElement() { };
	PatternClipboardElement(const CString &data, const CString &desc) : content(data), description(desc) { };
};


//====================
class PatternClipboard
//====================
{
	friend class PatternClipboardDialog;

public:

	enum PasteModes
	{
		pmOverwrite = 0,
		pmMixPaste,
		pmMixPasteIT,
		pmPasteFlood,
		pmPushForward,
	};

	// Clipboard index type
	typedef size_t clipindex_t;

protected:

	// The one and only pattern clipboard object.
	static PatternClipboard patternClipboard;

	// Active internal clipboard index
	clipindex_t activeClipboard;
	// Internal clipboard collection
	std::deque<PatternClipboardElement> clipboards;

public:

	PatternClipboard() : activeClipboard(0) { SetClipboardSize(1); };

	// Copy a range of patterns to both the system clipboard and the internal clipboard.
	bool Copy(CSoundFile &sndFile, ORDERINDEX first, ORDERINDEX last);
	// Copy a pattern selection to both the system clipboard and the internal clipboard.
	bool Copy(CSoundFile &sndFile, PATTERNINDEX pattern, PatternRect selection);
	// Try pasting a pattern selection from the system clipboard.
	bool Paste(CSoundFile &sndFile, ModCommandPos &pastePos, PasteModes mode, ORDERINDEX curOrder);
	// Try pasting a pattern selection from an internal clipboard.
	bool Paste(CSoundFile &sndFile, ModCommandPos &pastePos, PasteModes mode, ORDERINDEX curOrder, clipindex_t internalClipboard);
	// Copy one of the internal clipboards to the system clipboard.
	bool SelectClipboard(clipindex_t which);
	// Switch to the next internal clipboard.
	bool CycleForward();
	// Switch to the previous internal clipboard.
	bool CycleBackward();
	// Set the maximum number of internal clipboards.
	void SetClipboardSize(clipindex_t maxEntries);
	// Return the current number of clipboards.
	clipindex_t GetClipboardSize() const { return clipboards.size(); };

	// Get the pattern clipboard singleton
	static PatternClipboard &Instance() { return patternClipboard; }

protected:

	// Create the clipboard text for a pattern selection
	CString CreateClipboardString(CSoundFile &sndFile, PATTERNINDEX pattern, PatternRect selection);

	CString GetFileExtension(const char *ext) const;
	// Parse clipboard string and perform the pasting operation.
	bool HandlePaste(CSoundFile &sndFile, ModCommandPos &pastePos, PasteModes mode, const CString &data, ORDERINDEX curOrder);

	// System-specific clipboard functions
	bool ToSystemClipboard(const PatternClipboardElement &clipboard) { return ToSystemClipboard(clipboard.content); };
	bool ToSystemClipboard(const CString &data);
	bool FromSystemClipboard(CString &data);
};


//===========================================
class PatternClipboardDialog : public CDialog
//===========================================
{
protected:

	// The one and only pattern clipboard dialog object
	static PatternClipboardDialog patternClipboardDialog;

	PatternClipboard &clipboards;
	CSpinButtonCtrl numClipboardsSpin;
	CListBox clipList;
	int posX, posY;
	bool isLocked, isCreated;

public:

	PatternClipboardDialog(PatternClipboard &c);
	void UpdateList();
	void Show();
	void Toggle() { if(isCreated) OnCancel(); else Show(); }

	// Get the pattern clipboard dialog singleton
	static PatternClipboardDialog &Instance() { return patternClipboardDialog; }

protected:

	virtual void DoDataExchange(CDataExchange* pDX);
	DECLARE_MESSAGE_MAP();

	afx_msg void OnCancel();
	afx_msg void OnNumClipboardsChanged();
	afx_msg void OnSelectClipboard();
};
