import React, { useEffect, useState, useRef } from 'react';
import { GoogleGenerativeAI } from '@google/generative-ai';

const WAKE_WORD = 'he buddy';
const API_KEY = process.env.REACT_APP_GEMINI_API_KEY;

export default function GeminiChatbot() {
  const [isListening, setIsListening] = useState(false);
  const [status, setStatus] = useState('Initializing...');
  const [conversation, setConversation] = useState([]);
  const [isProcessing, setIsProcessing] = useState(false);
  const [error, setError] = useState(null);
  const [lastTranscript, setLastTranscript] = useState('');
  
  const recognitionRef = useRef(null);
  const genAIRef = useRef(null);
  const wakeWordTimeoutRef = useRef(null);
  
  useEffect(() => {
    console.log('API Key present:', !!API_KEY);
    console.log('API Key starts with AIzaSy:', API_KEY?.startsWith('AIzaSy'));
    
    if (!API_KEY) {
      setError('Missing REACT_APP_GEMINI_API_KEY in .env.local file');
      setStatus('Configuration Error');
      return;
    }
    
    try {
      genAIRef.current = new GoogleGenerativeAI(API_KEY);
      testApiKey();
    } catch (error) {
      setError('Failed to initialize Gemini AI: ' + error.message);
      setStatus('Initialization Error');
      return;
    }
    
    initSpeechRecognition();
  }, []);
  
  const testApiKey = async () => {
    try {
      const model = genAIRef.current.getGenerativeModel({ model: 'gemini-1.5-flash' });
      const result = await model.generateContent('Say hello in one word');
      const response = await result.response;
      const text = await response.text();
      console.log('API test successful:', text);
      setStatus('✅ Ready - Say "hey buddy" + your question');
      setError(null);
    } catch (error) {
      console.error('API test failed:', error);
      setError('API key test failed: ' + error.message);
      setStatus('❌ API Error - Check your key');
    }
  };
  
  const initSpeechRecognition = () => {
    if ('webkitSpeechRecognition' in window || 'SpeechRecognition' in window) {
      const SpeechRecognition = window.SpeechRecognition || window.webkitSpeechRecognition;
      recognitionRef.current = new SpeechRecognition();
      
      const recognition = recognitionRef.current;
      recognition.continuous = true;
      recognition.interimResults = true;
      recognition.lang = 'en-US';
      
      recognition.onstart = () => {
        setIsListening(true);
        setStatus('🎤 Listening for "hey buddy"...');
      };
      
      recognition.onresult = (event) => {
        let finalTranscript = '';
        let interimTranscript = '';
        
        for (let i = event.resultIndex; i < event.results.length; i++) {
          const transcript = event.results[i][0].transcript;
          if (event.results[i].isFinal) {
            finalTranscript += transcript;
          } else {
            interimTranscript += transcript;
          }
        }
        
        const fullTranscript = (finalTranscript + interimTranscript).toLowerCase();
        setLastTranscript(fullTranscript);
        
        if (fullTranscript.includes(WAKE_WORD)) {
          const wakeWordIndex = fullTranscript.lastIndexOf(WAKE_WORD);
          const commandPart = fullTranscript.substring(wakeWordIndex + WAKE_WORD.length).trim();
          
          console.log('Wake word detected! Command:', commandPart);
          setStatus('🔥 Wake word detected! Processing...');
          
          if (wakeWordTimeoutRef.current) {
            clearTimeout(wakeWordTimeoutRef.current);
          }
          
          if (commandPart.length > 2) {
            handleVoiceCommand(commandPart);
          } else {
            wakeWordTimeoutRef.current = setTimeout(() => {
              const latestTranscript = (finalTranscript + interimTranscript).toLowerCase();
              const latestCommand = latestTranscript.substring(latestTranscript.lastIndexOf(WAKE_WORD) + WAKE_WORD.length).trim();
              
              if (latestCommand.length > 2) {
                handleVoiceCommand(latestCommand);
              } else {
                setStatus('❌ No command after wake word. Try again.');
                setTimeout(() => setStatus('🎤 Listening for "hey buddy"...'), 2000);
              }
            }, 2500);
          }
        }
      };
      
      recognition.onerror = (event) => {
        console.error('Speech error:', event.error);
        setStatus('❌ Speech error: ' + event.error);
        setTimeout(() => {
          if (!error) startListening();
        }, 2000);
      };
      
      recognition.onend = () => {
        setIsListening(false);
        if (!error && !status.includes('Error')) {
          setTimeout(startListening, 1000);
        }
      };
      
      if (!error) {
        startListening();
      }
    } else {
      setStatus('❌ Speech recognition not supported');
    }
  };
  
  const startListening = () => {
    if (recognitionRef.current && !isListening && !error) {
      try {
        recognitionRef.current.start();
      } catch (error) {
        console.log('Recognition already running');
      }
    }
  };
  
  const stopListening = () => {
    if (recognitionRef.current) {
      recognitionRef.current.stop();
      setIsListening(false);
      setStatus('Stopped');
    }
    if (wakeWordTimeoutRef.current) {
      clearTimeout(wakeWordTimeoutRef.current);
    }
  };
  
  const handleVoiceCommand = async (command) => {
    console.log('Processing command:', command);
    
    if (!genAIRef.current) {
      speak('AI system is not ready');
      return;
    }
    
    if (wakeWordTimeoutRef.current) {
      clearTimeout(wakeWordTimeoutRef.current);
    }
    
    setIsProcessing(true);
    setStatus('🤖 Thinking...');
    
    const cleanCommand = command.replace(/[^\w\s]/g, '').trim();
    const userMessage = { type: 'user', text: cleanCommand };
    setConversation(prev => [...prev, userMessage]);
    
    try {
      const model = genAIRef.current.getGenerativeModel({ 
        model: 'gemini-1.5-flash',
        generationConfig: {
          maxOutputTokens: 150,
          temperature: 0.7,
        }
      });
      
      const prompt = `You are buddy, a helpful AI tutor. Answer briefly and conversationally: ${cleanCommand}`;
      
      const result = await model.generateContent(prompt);
      const response = await result.response;
      const text = response.text();
      
      console.log('Response:', text);
      
      const aiMessage = { type: 'ai', text };
      setConversation(prev => [...prev, aiMessage]);
      
      speak(text);
      setStatus('✅ Answered! Say "hey buddy" for more');
      
    } catch (error) {
      console.error('API error:', error);
      const errorMsg = 'Sorry, I had trouble with that.';
      setConversation(prev => [...prev, { type: 'ai', text: errorMsg }]);
      speak(errorMsg);
      setStatus('❌ Error. Try again.');
    }
    
    setIsProcessing(false);
    setTimeout(() => {
      setStatus('🎤 Listening for "hey buddy"...');
      if (!isListening) startListening();
    }, 3000);
  };
  
  const speak = (text) => {
    if ('speechSynthesis' in window) {
      speechSynthesis.cancel();
      const utterance = new SpeechSynthesisUtterance(text);
      utterance.rate = 0.8;
      utterance.pitch = 1.0;
      speechSynthesis.speak(utterance);
    }
  };
  
  const testCommand = () => {
    handleVoiceCommand('what is photosynthesis');
  };
  
  if (error) {
    return (
      <div className="gemini-chatbot">
        <div className="error-state">
          <div className="error-icon">⚠️</div>
          <div className="error-message">{error}</div>
          <button onClick={() => window.location.reload()}>Retry</button>
        </div>
      </div>
    );
  }
  
  return (
    <div className="gemini-chatbot">
      <div className="chatbot-header">
        <div className="status">
          <span className={`status-indicator ${isListening ? 'listening' : ''} ${isProcessing ? 'processing' : ''}`}>
            {isProcessing ? '🤖' : isListening ? '🎤' : '🔇'}
          </span>
          <span className="status-text">{status}</span>
        </div>
        <div className="controls">
          <button onClick={startListening} disabled={isListening}>Start</button>
          <button onClick={stopListening} disabled={!isListening}>Stop</button>
          <button onClick={() => setConversation([])}>Clear</button>
          <button onClick={testCommand} className="test-btn">Test</button>
        </div>
      </div>
      
      <div className="wake-word-info">
        <div className="small">Say <strong>"hey buddy"</strong> + question</div>
        <div className="small">Example: "Hey buddy, what is gravity?"</div>
        {lastTranscript && (
          <div className="small transcript">Heard: "{lastTranscript.slice(-50)}"</div>
        )}
      </div>
      
      <div className="conversation">
        {conversation.length === 0 && (
          <div className="empty-state">
            <div className="small">Try: "Hey buddy, explain photosynthesis"</div>
          </div>
        )}
        {conversation.map((message, index) => (
          <div key={index} className={`message ${message.type}`}>
            <div className="message-label">{message.type === 'user' ? 'You:' : 'buddy:'}</div>
            <div className="message-text">{message.text}</div>
          </div>
        ))}
        {isProcessing && (
          <div className="message ai">
            <div className="message-label">buddy:</div>
            <div className="message-text typing">Thinking...</div>
          </div>
        )}
      </div>
    </div>
  );
}
