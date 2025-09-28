import React, { useEffect, useState } from 'react';
import * as Vosk from 'vosk-browser';

export default function VoiceCommands() {
  const [status, setStatus] = useState('idle');
  const [lastText, setLastText] = useState('');

  useEffect(() => {
    let audioContext, recognizerNode, source, recognizer, stopped = false;
    (async () => {
      try {
        setStatus('loading');
        const model = await Vosk.createModel('/model.tar.gz');
        recognizer = new model.KaldiRecognizer();
        recognizer.on('result', (m) => setLastText(m.result.text || ''));
        recognizer.on('partialresult', (m) => setLastText(m.result.partial || ''));
        const media = await navigator.mediaDevices.getUserMedia({
          audio: { echoCancellation: true, noiseSuppression: true, channelCount: 1, sampleRate: 16000 }
        });
        audioContext = new (window.AudioContext || window.webkitAudioContext)();
        recognizerNode = audioContext.createScriptProcessor(4096, 1, 1);
        recognizerNode.onaudioprocess = (e) => {
          if (!stopped) {
            try { recognizer.acceptWaveform(e.inputBuffer); } catch (err) { /* ignore */ }
          }
        };
        source = audioContext.createMediaStreamSource(media);
        source.connect(recognizerNode);
        setStatus('listening');
      } catch (e) {
        console.error(e);
        setStatus('error');
      }
    })();
    return () => {
      stopped = true;
      try { source && source.disconnect(); } catch {}
      try { recognizerNode && recognizerNode.disconnect(); } catch {}
      try { audioContext && audioContext.close(); } catch {}
    };
  }, []);

  return (
    <div>
      <div className="small">Status: {status}</div>
      <div className="small">Heard: {lastText}</div>
    </div>
  );
}
